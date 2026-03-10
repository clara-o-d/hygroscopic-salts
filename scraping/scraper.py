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
#  CONFIG
# ─────────────────────────────────────────────

TEST_MODE  = False  # True → small sample only; False → full run
N_WORKERS  = 10     # parallel browser instances (ignored in TEST_MODE)

URL = "https://members.cere.dk/software/edatabase/"

# ─────────────────────────────────────────────
#  ION DEFINITIONS  (symbol → charge magnitude)
# ─────────────────────────────────────────────

CATIONS = {
    "Li":  1,
    "Na":  1,
    "K":   1,
    "Cs":  1,
    "NH4": 1,
    "Mg":  2,
    "Ca":  2,
    "Sr":  2,
    "Ba":  2,
}

ANIONS = {
    "F":    1,
    "Cl":   1,
    "Br":   1,
    "I":    1,
    "OH":   1,
    "NO3":  1,
    "HCO3": 1,
    "SO4":  2,
    "CO3":  2,
}

POLYATOMIC = {"NH4", "NO3", "OH", "HCO3", "SO4", "CO3", "PO4"}

# ─────────────────────────────────────────────
#  FORMULA GENERATION
# ─────────────────────────────────────────────

def make_formula(cation: str, cat_charge: int, anion: str, an_charge: int) -> str:
    lcm   = (cat_charge * an_charge) // gcd(cat_charge, an_charge)
    n_cat = lcm // cat_charge
    n_an  = lcm // an_charge

    cat_str = (f"({cation}){n_cat}" if cation in POLYATOMIC and n_cat > 1
               else cation + (str(n_cat) if n_cat > 1 else ""))
    an_str  = (f"({anion}){n_an}"  if anion  in POLYATOMIC and n_an  > 1
               else anion  + (str(n_an)  if n_an  > 1 else ""))
    return cat_str + an_str


single_salts   = [make_formula(c, cc, a, ac)
                  for c, cc in CATIONS.items()
                  for a, ac in ANIONS.items()]
two_comp_pairs = list(combinations(single_salts, 2))

# ─────────────────────────────────────────────
#  TEST SUBSET
# ─────────────────────────────────────────────

TEST_SINGLE = ["NaCl", "CaCl2", "KCl", "MgCl2"]
TEST_TWO    = [("NaCl", "KCl"), ("MgCl2", "CaCl2"), ("NaCl", "CaCl2"), ("KCl", "MgCl2")]

if TEST_MODE:
    single_salts   = TEST_SINGLE
    two_comp_pairs = TEST_TWO

n_workers = 1 if TEST_MODE else N_WORKERS
print(f"{'[TEST]' if TEST_MODE else '[FULL]'} "
      f"{len(single_salts)} singles, {len(two_comp_pairs)} pairs "
      f"across {n_workers} workers")

# ─────────────────────────────────────────────
#  WORK DISTRIBUTION  (interleaved so each worker gets a diverse slice)
# ─────────────────────────────────────────────

def interleaved_chunks(lst, n):
    """Split lst into n chunks by interleaving (worker i gets items i, i+n, i+2n, …)."""
    return [lst[i::n] for i in range(n)]

single_chunks = interleaved_chunks(single_salts,   n_workers)
pair_chunks   = interleaved_chunks(two_comp_pairs,  n_workers)

# ─────────────────────────────────────────────
#  SHARED STATE
# ─────────────────────────────────────────────

all_rows  = []
rows_lock = Lock()
bar_lock  = Lock()
stem      = "cere_vle_test" if TEST_MODE else "cere_vle_data"

s_counts = {"data": 0, "no_data": 0, "error": 0}
p_counts = {"data": 0, "no_data": 0, "error": 0}

single_bar = tqdm(total=len(single_salts),   desc="Singles", unit="salt",
                  position=0, leave=True, dynamic_ncols=True, disable=TEST_MODE)
pair_bar   = tqdm(total=len(two_comp_pairs), desc="Pairs  ", unit="pair",
                  position=1, leave=True, dynamic_ncols=True, disable=TEST_MODE)

def _save_unlocked():
    if not all_rows:
        return
    df = pd.DataFrame(all_rows)
    df.to_csv(f"{stem}.csv",  index=False)
    df.to_json(f"{stem}.json", orient="records", indent=2)

def save_progress(msg=""):
    with rows_lock:
        _save_unlocked()
        tqdm.write(f"  💾 Saved {len(all_rows)} rows{(' — ' + msg) if msg else ''}")

# ─────────────────────────────────────────────
#  WORKER CLASS
# ─────────────────────────────────────────────

class Worker:
    def __init__(self, wid: int, singles: list, pairs: list):
        self.wid     = wid
        self.singles = singles
        self.pairs   = pairs
        self.driver  = None
        self.wait    = None
        self._launch_browser()

    # ── Browser lifecycle ──────────────────────

    def _launch_browser(self):
        self.driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()))
        self.wait   = WebDriverWait(self.driver, 15)

    def _restart_browser(self):
        tqdm.write(f"  [W{self.wid}] ⚠ Browser crashed — saving & restarting...")
        save_progress(f"worker {self.wid} crash")
        try:
            self.driver.quit()
        except Exception:
            pass
        time.sleep(3)
        self._launch_browser()
        for delay in [5, 15, 30, 60]:
            try:
                self.driver.get(URL)
                time.sleep(2)
                return
            except Exception as e:
                tqdm.write(f"  [W{self.wid}] ⚠ Net error ({e.__class__.__name__}) — retry in {delay}s...")
                time.sleep(delay)
        self.driver.get(URL)
        time.sleep(2)

    @staticmethod
    def _is_dead(exc):
        msg = str(exc).lower()
        return "invalid session id" in msg or "no such window" in msg

    # ── Page interaction ───────────────────────

    def _select_prop_and_ncomp(self, n_comp: int):
        prop = Select(self.wait.until(
            EC.presence_of_element_located((By.TAG_NAME, "select"))
        ))
        prop.select_by_value("vle")
        time.sleep(1)
        selects = self.driver.find_elements(By.TAG_NAME, "select")
        Select(selects[1]).select_by_value(str(n_comp))
        time.sleep(0.5)

    def _click_search(self):
        self.driver.find_element(
            By.XPATH, "//button[contains(text(),'Search')]"
        ).click()

    def _collect(self, label: str) -> int:
        local = []
        for table in self.driver.find_elements(By.TAG_NAME, "table"):
            rows = table.find_elements(By.TAG_NAME, "tr")
            if len(rows) < 2:
                continue
            headers = [c.text for c in rows[0].find_elements(By.TAG_NAME, "th")]
            if not headers:
                continue
            for r in rows[1:]:
                cols = [c.text for c in r.find_elements(By.TAG_NAME, "td")]
                if len(cols) == len(headers):
                    row_dict = dict(zip(headers, cols))
                    row_dict["salt"] = label
                    local.append(row_dict)
        with rows_lock:
            all_rows.extend(local)
        return len(local)

    # ── Core scrape (1 or 2 components) ───────

    def _scrape(self, n_comp: int, formulas: tuple, label: str) -> int:
        for attempt in range(2):
            try:
                self.driver.get(URL)
                time.sleep(1)
                self._select_prop_and_ncomp(n_comp)
                inputs = self.driver.find_elements(By.XPATH, "//input[@type='text']")
                for i, f in enumerate(formulas):
                    inputs[i].clear()
                    inputs[i].send_keys(f)
                self._click_search()
                time.sleep(15)
                return self._collect(label)
            except Exception as e:
                if attempt == 0 and self._is_dead(e):
                    self._restart_browser()
                else:
                    tqdm.write(f"  [W{self.wid}] ✗ {label} → {e.__class__.__name__}: {str(e)[:80]}")
                    return -1   # signal error
        return -1

    # ── Main run loop ──────────────────────────

    def run(self):
        self.driver.get(URL)
        time.sleep(2)

        # Singles
        for salt in self.singles:
            added = self._scrape(1, (salt,), salt)
            with bar_lock:
                if added > 0:
                    s_counts["data"]    += 1
                elif added == 0:
                    s_counts["no_data"] += 1
                else:
                    s_counts["error"]   += 1
                single_bar.update(1)
                single_bar.set_postfix(s_counts)
            if TEST_MODE:
                tqdm.write(f"  [W{self.wid}] ✓ {salt}  (+{added})")

        # Pairs
        for salt1, salt2 in self.pairs:
            label = f"{salt1}+{salt2}"
            added = self._scrape(2, (salt1, salt2), label)
            with bar_lock:
                if added > 0:
                    p_counts["data"]    += 1
                elif added == 0:
                    p_counts["no_data"] += 1
                else:
                    p_counts["error"]   += 1
                pair_bar.update(1)
                pair_bar.set_postfix(p_counts)
            if TEST_MODE:
                tqdm.write(f"  [W{self.wid}] ✓ {label}  (+{added})")

        try:
            self.driver.quit()
        except Exception:
            pass

# ─────────────────────────────────────────────
#  LAUNCH WORKERS
# ─────────────────────────────────────────────

workers = [Worker(i, single_chunks[i], pair_chunks[i]) for i in range(n_workers)]
threads = [Thread(target=w.run, daemon=True) for w in workers]

for t in threads:
    t.start()
for t in threads:
    t.join()

single_bar.close()
pair_bar.close()

# ─────────────────────────────────────────────
#  SAVE FINAL
# ─────────────────────────────────────────────

df = pd.DataFrame(all_rows)
df.to_csv(f"{stem}.csv",  index=False)
df.to_json(f"{stem}.json", orient="records", indent=2)
print(f"\nDone. Saved {len(df)} rows → {stem}.csv / {stem}.json")
print(f"  Singles → {s_counts}")
print(f"  Pairs   → {p_counts}")
