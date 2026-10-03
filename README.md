# 🦋 NEO-CRISPR Core Engine

**A modern, fully local CRISPR guide design tool built for biologists working with non-model organisms.**

NEO-CRISPR is a standalone desktop application that runs entirely on your own computer. You do not need to upload your sensitive sequences to a web server, and there are no cloud computing costs. It features a modern interface with smart guide ranking, interactive sequence highlighting, and rapid genome-wide off-target scanning.

While designed with a focus on *Lepidoptera*, this tool works perfectly for **any organism**.

<img width="2264" height="1324" alt="neo_crispr" src="https://github.com/user-attachments/assets/321fe5c9-e76b-47fc-bc5a-b1e95c7437da" />


---

## ⚡ Features

* **Absolute Privacy:** Everything runs locally on your machine. The interface ships
  with the app (no CDN, no fonts fetched from the internet), so it also works fully offline.
* **Both Strands:** Candidate guides are found on **both** the forward and reverse strands
  of your target sequence.
* **Genome-Wide Off-Target Scan:** Every scaffold of the genome is scanned on **both strands**
  for off-targets (up to 3 mismatches), tolerating both the canonical `NGG` PAM and the weaker
  `NAG` PAM. Each guide's own on-target site is excluded, and the PAM-proximal "seed" mismatches
  are reported separately because they matter most for specificity.
* **Gene Annotation:** If you provide a `.gtf` file, each off-target is labelled with the gene
  it lands in.
* **Plug-and-Play Genomes:** No coding required to add a new species. Just drop your `.fa`
  (FASTA) file into a folder (a `.gtf` is optional, for gene annotation), and the app detects it
  automatically.
* **Interactive Sequence Map:** Hover over any generated guide to instantly see exactly where it
  binds on your target sequence, with the spacer and PAM highlighted separately.
* **Hardware Accelerated:** Uses tensor math (PyTorch) for the off-target scan, using your
  computer's GPU if available, or seamlessly running on the CPU. The genome is encoded once and
  cached, and every candidate guide is scanned in a single batched pass.

> **A note on scoring:** the per-guide **score** is a transparent rule of thumb (it rewards
> mid-range GC content and a 3′ G). It is *not* a validated on-target efficiency predictor such
> as Rule Set 2 or DeepSpCas9 — use it as a tie-breaker, not as ground truth.

---

## 🛠️ Getting Started 

You do not need to be a computer scientist to run this tool. Just follow these steps once to set up your computer.

### Step 1: Download the App
1. At the top of this GitHub page, click the green **"<> Code"** button.
2. Select **"Download ZIP"**.
3. Extract the downloaded ZIP file to a permanent location on your computer (for example, your Documents folder or Desktop).

### Step 2: Install Python
If you do not already have Python installed on your computer:
1. Go to [Python.org](https://www.python.org/downloads/) and download the latest version for your operating system (Windows or Mac).
2. Run the installer. **CRITICAL:** When the installer opens, make sure to check the box at the bottom that says **"Add python.exe to PATH"** before clicking Install.


### Step 3: Install Required Libraries
This app relies on a few standard scientific libraries to do the heavy lifting.

1. Open your computer's terminal:
   * **Windows:** Press the Start button, type `cmd` or `PowerShell`, and hit Enter.
   * **Mac:** Open Spotlight (Command + Space), type `Terminal`, and hit Enter.
2. Copy and paste this exact command into the terminal and press Enter:
   ```bash
   pip install -r requirements.txt
   ```
   (Or install the packages directly: `pip install fastapi uvicorn pyfaidx torch pandas`.)

### Step 4: Set Up Your Organism Data
Move to the folder where you have unzipped the app. Inside the `data/` folder, create a folder
with your organism's scientific name (for example, `bicyclus_anynana`) and add your genome FASTA
as `genome.fa`. A GTF annotation (`annotation.gtf`) is **optional** — add it if you want each
off-target labelled with the gene it falls in.

### Step 5: Start the App
3. Copy and paste this exact command into the terminal (in the same location where app.py is present) and press Enter:
   ```bash
   python app.py
   ```
   <img width="1893" height="889" alt="Screenshot 2026-07-17 075109" src="https://github.com/user-attachments/assets/d6385edc-5dad-4481-962c-c7372cce269a" />

---

## 🧪 For Developers

The scientific core lives in `engine.py` (sequence cleaning, both-strand guide discovery,
genome-wide off-target scanning, GTF annotation) and is independent of the web server in
`app.py`. The browser UI is served from `static/`.

Run the test suite with:

```bash
pip install -r requirements.txt
pytest
```

---

Designed by [Tirtha Das Banerjee](https://tirthadasbanerjee.com/).
