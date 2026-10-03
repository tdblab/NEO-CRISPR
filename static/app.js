"use strict";

let speciesMeta = {};

function esc(s) {
  return String(s).replace(/[&<>"']/g, (c) => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
  ));
}

document.addEventListener("DOMContentLoaded", async () => {
  const select = document.getElementById("speciesSelect");
  try {
    const res = await fetch("/api/species");
    const species = await res.json();
    if (!species.length) {
      select.innerHTML = '<option value="">No genomes found in data/ folder</option>';
      return;
    }
    species.forEach((s) => { speciesMeta[s.id] = s; });
    select.innerHTML = species
      .map((s) => `<option value="${esc(s.id)}">${esc(s.name)}</option>`)
      .join("");
    updateAnnotationHint();
    select.addEventListener("change", updateAnnotationHint);
    document.getElementById("seqInput").value =
      "ATGCCGCGCGTCGTGCCCGACCAGAGAAGCAAGTTCGAGAACGG";
  } catch (e) {
    select.innerHTML = '<option value="">Could not reach local engine</option>';
  }
});

function updateAnnotationHint() {
  const id = document.getElementById("speciesSelect").value;
  const meta = speciesMeta[id];
  const hint = document.getElementById("annotationHint");
  if (!meta) { hint.textContent = ""; return; }
  hint.textContent = meta.annotated
    ? "GTF annotation detected: off-targets will be labelled with their gene."
    : "No GTF found: off-targets will not be gene-annotated (add annotation.gtf to enable).";
}

function renderSequenceMap(seq) {
  const container = document.getElementById("sequenceMap");
  container.innerHTML = Array.from(seq)
    .map((ch, i) => `<span class="nt" id="nt-${i}">${ch}</span>`)
    .join("");
  document.getElementById("sequenceViewer").classList.remove("hidden");
}

window.highlightGuide = function (start, end, strand) {
  // PAM is 3 nt; its forward-coord position depends on strand.
  const pamStart = strand === "+" ? end - 3 : start;
  const pamEnd = strand === "+" ? end : start + 3;
  for (let i = start; i < end; i++) {
    const el = document.getElementById(`nt-${i}`);
    if (!el) continue;
    if (i >= pamStart && i < pamEnd) el.classList.add("hl-pam");
    else el.classList.add("hl");
  }
};

window.clearHighlight = function (start, end) {
  for (let i = start; i < end; i++) {
    const el = document.getElementById(`nt-${i}`);
    if (el) el.classList.remove("hl", "hl-pam");
  }
};

async function processSequence() {
  const sequence = document.getElementById("seqInput").value;
  const speciesId = document.getElementById("speciesSelect").value;
  if (!sequence.trim() || !speciesId) return;

  document.getElementById("emptyState").classList.add("hidden");
  const results = document.getElementById("resultsContent");
  results.classList.remove("hidden");
  results.innerHTML =
    '<div class="status loading">Encoding genome &amp; scanning for off-targets (whole genome, both strands)...</div>';
  document.getElementById("sequenceViewer").classList.add("hidden");

  let data;
  try {
    const response = await fetch("/api/design", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sequence, species_id: speciesId }),
    });
    data = await response.json();
  } catch (e) {
    results.innerHTML = '<div class="status warn">Request failed. Is the local engine running?</div>';
    return;
  }

  if (data.error) {
    results.innerHTML = `<div class="status warn">${esc(data.error)}</div>`;
    return;
  }
  if (data.device) {
    document.getElementById("engineState").textContent =
      "LOCAL " + data.device + " READY";
  }

  if (!data.candidates || !data.candidates.length) {
    results.innerHTML =
      '<div class="status warn">No SpCas9 (NGG) target sites found on either strand of this sequence.</div>';
    return;
  }

  renderSequenceMap(data.input_sequence);
  results.innerHTML = data.candidates.map((c) => renderCard(c, data.annotated)).join("");
}

function renderCard(c, annotated) {
  const off = c.off_targets;
  let offHtml;
  if (!off.length) {
    offHtml = '<div class="clean">&#10003; No off-targets found genome-wide (&le;3 mismatches)</div>';
  } else {
    const geneCol = annotated ? "<th>Gene</th>" : "";
    const rows = off.map((ot) => `
      <tr>
        <td>${esc(ot.loc)}</td>
        <td>${esc(ot.seq)}</td>
        <td>${esc(ot.pam)}</td>
        <td class="right">${ot.seed_mm}</td>
        <td class="mm">${ot.mm}</td>
        ${annotated ? `<td class="gene">${ot.gene ? esc(ot.gene) : "&mdash;"}</td>` : ""}
      </tr>`).join("");
    offHtml = `
      <table>
        <thead><tr>
          <th>Locus (strand)</th><th>Alignment</th><th>PAM</th>
          <th class="right">Seed MM</th><th class="right">Total MM</th>${geneCol}
        </tr></thead>
        <tbody>${rows}</tbody>
      </table>`;
  }

  return `
    <div class="card ${c.is_best ? "best" : ""}"
         onmouseenter="highlightGuide(${c.start}, ${c.end}, '${c.strand}')"
         onmouseleave="clearHighlight(${c.start}, ${c.end})">
      ${c.is_best ? '<div class="badge">&#127942; Recommended Guide</div>' : ""}
      <div class="card-head">
        <span class="seq mono">${esc(c.spacer)}<span class="pam">${esc(c.pam)}</span></span>
        <div class="tags">
          <span class="tag strand">Strand ${esc(c.strand)}</span>
          <span class="tag">GC ${esc(c.gc)}</span>
          <span class="tag eff">Score ${esc(c.score)}</span>
        </div>
      </div>
      <div class="offsec">
        <div class="title">Genome-wide off-target profile (${off.length} site${off.length === 1 ? "" : "s"}):</div>
        ${offHtml}
      </div>
    </div>`;
}
