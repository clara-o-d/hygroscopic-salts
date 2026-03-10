from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

from itertools import combinations
from math import gcd
from threading import Thread, Lock

import pandas as pd
import time
from tqdm import tqdm


# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────

TEST_MODE = False
N_WORKERS = 20

URL = "https://members.cere.dk/software/edatabase/"

MISC_PROPS = [
    ("water-activity", "water_activity"),
    ("osmotic-coefficients", "osmotic_coeff"),
    ("activity-coefficients", "activity_coeff"),
]


# ─────────────────────────────────────────────
# IONS
# ─────────────────────────────────────────────

CATIONS = {
    "Li": 1,
    "Na": 1,
    "K": 1,
    "Cs": 1,
    "NH4": 1,
    "Mg": 2,
    "Ca": 2,
    "Sr": 2,
    "Ba": 2,
}

ANIONS = {
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "OH": 1,
    "NO3": 1,
    "HCO3": 1,
    "SO4": 2,
    "CO3": 2,
}

POLYATOMIC = {"NH4", "NO3", "OH", "HCO3", "SO4", "CO3"}


# ─────────────────────────────────────────────
# FORMULA GENERATION
# ─────────────────────────────────────────────

def make_formula(cation, cat_charge, anion, an_charge):

    lcm = (cat_charge * an_charge) // gcd(cat_charge, an_charge)

    n_cat = lcm // cat_charge
    n_an = lcm // an_charge

    cat_str = (
        f"({cation}){n_cat}"
        if cation in POLYATOMIC and n_cat > 1
        else cation + (str(n_cat) if n_cat > 1 else "")
    )

    an_str = (
        f"({anion}){n_an}"
        if anion in POLYATOMIC and n_an > 1
        else anion + (str(n_an) if n_an > 1 else "")
    )

    return cat_str + an_str


single_salts = [
    make_formula(c, cc, a, ac)
    for c, cc in CATIONS.items()
    for a, ac in ANIONS.items()
]

two_comp_pairs = list(combinations(single_salts, 2))


TEST_SINGLE = ["NaCl", "CaCl2", "KCl", "MgCl2"]

TEST_TWO = [
    ("NaCl", "KCl"),
    ("MgCl2", "CaCl2"),
    ("NaCl", "CaCl2"),
    ("KCl", "MgCl2"),
]


if TEST_MODE:

    single_salts = TEST_SINGLE
    two_comp_pairs = TEST_TWO


n_workers = 1 if TEST_MODE else N_WORKERS

print(
    f"{'[TEST]' if TEST_MODE else '[FULL]'} "
    f"{len(single_salts)} singles, "
    f"{len(two_comp_pairs)} pairs "
    f"across {n_workers} workers"
)


# ─────────────────────────────────────────────
# SPLIT WORK
# ─────────────────────────────────────────────

def interleaved_chunks(lst, n):
    return [lst[i::n] for i in range(n)]


single_chunks = interleaved_chunks(single_salts, n_workers)
pair_chunks = interleaved_chunks(two_comp_pairs, n_workers)


# ─────────────────────────────────────────────
# SHARED STATE
# ─────────────────────────────────────────────

all_rows = []

rows_lock = Lock()
bar_lock = Lock()

stem = "cere_misc_test" if TEST_MODE else "cere_misc_data"


s_counts = {"data": 0, "no_data": 0, "error": 0}
p_counts = {"data": 0, "no_data": 0, "error": 0}


single_bar = tqdm(
    total=len(single_salts),
    desc="Singles",
    unit="salt",
    position=0,
    leave=True,
)

pair_bar = tqdm(
    total=len(two_comp_pairs),
    desc="Pairs",
    unit="pair",
    position=1,
    leave=True,
)


# ─────────────────────────────────────────────
# WORKER
# ─────────────────────────────────────────────


class Worker:

    def __init__(self, wid, singles, pairs):

        self.wid = wid
        self.singles = singles
        self.pairs = pairs

        self.driver = None
        self.wait = None

        self.launch()


    def launch(self):

        self.driver = webdriver.Chrome(
            service=Service(
                ChromeDriverManager().install()
            )
        )

        self.wait = WebDriverWait(
            self.driver,
            15,
        )


    def restart(self):

        try:
            self.driver.quit()
        except:
            pass

        time.sleep(2)

        self.launch()


    def select_misc(self, n_comp, prop_value):

        selects = self.wait.until(
            EC.presence_of_all_elements_located(
                (By.TAG_NAME, "select")
            )
        )

        Select(selects[0]).select_by_value("misc")

        time.sleep(1)

        selects = self.driver.find_elements(
            By.TAG_NAME,
            "select",
        )

        Select(selects[1]).select_by_value(
            prop_value
        )

        Select(selects[2]).select_by_value(
            str(n_comp)
        )


    def click_search(self):

        self.driver.find_element(
            By.XPATH,
            "//button[contains(text(),'Search')]"
        ).click()


    def collect(self, label, prop):

        local = []

        tables = self.driver.find_elements(
            By.TAG_NAME,
            "table"
        )

        for table in tables:

            rows = table.find_elements(
                By.TAG_NAME,
                "tr"
            )

            if len(rows) < 2:
                continue

            headers = [
                c.text
                for c in rows[0].find_elements(
                    By.TAG_NAME,
                    "th"
                )
            ]

            if not headers:
                continue

            for r in rows[1:]:

                cols = [
                    c.text
                    for c in r.find_elements(
                        By.TAG_NAME,
                        "td"
                    )
                ]

                if len(cols) == len(headers):

                    d = dict(
                        zip(headers, cols)
                    )

                    d["salt"] = label
                    d["property_type"] = prop

                    local.append(d)

        with rows_lock:
            all_rows.extend(local)

        return len(local)


    def scrape(self, n_comp, formulas, label):

        total = 0

        for prop_value, prop_name in MISC_PROPS:

            for attempt in range(2):

                try:

                    self.driver.get(URL)

                    time.sleep(1)

                    self.select_misc(
                        n_comp,
                        prop_value,
                    )

                    inputs = self.driver.find_elements(
                        By.XPATH,
                        "//input[@type='text']"
                    )

                    for i, f in enumerate(
                        formulas
                    ):
                        inputs[i].clear()
                        inputs[i].send_keys(f)

                    self.click_search()

                    time.sleep(10)

                    added = self.collect(
                        label,
                        prop_name,
                    )

                    total += added

                    break

                except Exception as e:

                    if attempt == 0:
                        self.restart()
                    else:
                        return -1

        return total


    def run(self):

        self.driver.get(URL)

        time.sleep(2)


        for salt in self.singles:

            added = self.scrape(
                1,
                (salt,),
                salt,
            )

            with bar_lock:

                if added > 0:
                    s_counts["data"] += 1
                elif added == 0:
                    s_counts["no_data"] += 1
                else:
                    s_counts["error"] += 1

                single_bar.update(1)


        for s1, s2 in self.pairs:

            label = f"{s1}+{s2}"

            added = self.scrape(
                2,
                (s1, s2),
                label,
            )

            with bar_lock:

                if added > 0:
                    p_counts["data"] += 1
                elif added == 0:
                    p_counts["no_data"] += 1
                else:
                    p_counts["error"] += 1

                pair_bar.update(1)


        try:
            self.driver.quit()
        except:
            pass


# ─────────────────────────────────────────────
# RUN
# ─────────────────────────────────────────────


workers = [
    Worker(
        i,
        single_chunks[i],
        pair_chunks[i],
    )
    for i in range(n_workers)
]


threads = [
    Thread(
        target=w.run,
        daemon=True,
    )
    for w in workers
]


for t in threads:
    t.start()

for t in threads:
    t.join()


single_bar.close()
pair_bar.close()


df = pd.DataFrame(all_rows)

df.to_csv(
    f"{stem}.csv",
    index=False,
)

df.to_json(
    f"{stem}.json",
    orient="records",
    indent=2,
)


print(
    f"\nDone. Saved {len(df)} rows"
)