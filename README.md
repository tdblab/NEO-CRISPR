# 🦋 NEO-CRISPR Core Engine

**A modern, fully local CRISPR guide design tool built for biologists working with non-model organisms.**

NEO-CRISPR is a standalone desktop application that runs entirely on your own computer. You do not need to upload your sensitive sequences to a web server, and there are no cloud computing costs. It features a modern interface with smart guide ranking, interactive sequence highlighting, and rapid genome-wide off-target scanning.

While designed with a focus on *Lepidoptera*, this tool works perfectly for **any organism**.

---

## ⚡ Features

* **Absolute Privacy:** Everything runs locally on your machine.
* **Plug-and-Play Genomes:** No coding required to add a new species. Just drop your `.fa` (FASTA) and `.gtf` annotation files into a folder, and the app will automatically detect and index them.
* **Smart Ranking:** Automatically ranks candidate guides to give you the best options first (prioritizing 0 off-targets and highest efficiency).
* **Interactive Sequence Map:** Hover over any generated guide to instantly see exactly where it binds on your target sequence.
* **Hardware Accelerated:** Uses advanced tensor math (PyTorch) to perform rapid off-target scanning, utilizing your computer's GPU if available, or seamlessly running on the CPU.

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
### Step 3: Install Required Libraries
This app relies on a few standard scientific libraries to do the heavy lifting.

1. Open your computer's terminal:
   * **Windows:** Press the Start button, type `cmd` or `PowerShell`, and hit Enter.
   * **Mac:** Open Spotlight (Command + Space), type `Terminal`, and hit Enter.
2. Copy and paste this exact command into the terminal and press Enter:
   ```bash
   pip install fastapi uvicorn pyfaidx torch pandas
   ```

### Step 4: Set Up Your Organism Data
Move to the folder where you have unzipped the app. Create a folder with your organism's scientific name (for example, `bicyclus_anynana`). Add the `genome.fa` and `annotation.gtf` for the organism inside that folder.

### Step 5: Start the App
3. Copy and paste this exact command into the terminal and press Enter:
   ```bash
   python app.py
   ```
Designed by [Tirtha Das Banerjee](https://tirthadasbanerjee.com/).
