/**
 * Typeahead suggestions for the single-search field (index.html).
 *
 * As the user types (2+ characters), fetches the best-matching proteins from
 * the /suggest endpoint and lists them under the input as
 * "GENE  Protein name  ACCESSION". Picking one (click, or arrow keys + Enter)
 * fills in its accession and submits the search, i.e. an exact lookup.
 * Pressing Enter without picking anything submits the typed text as before.
 */
(function () {
  'use strict';

  const input = document.getElementById('query');
  const list = document.getElementById('query-suggestions');
  if (!input || !list || !input.dataset.suggestUrl) return;

  const MIN_CHARS = 2;
  const DEBOUNCE_MS = 150;

  let items = [];          // current suggestions
  let active = -1;         // index of the keyboard-highlighted item
  let timer = null;
  let controller = null;   // aborts a superseded fetch

  function escapeHtml(s) {
    return String(s == null ? '' : s).replace(/[&<>"']/g, c => ({
      '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
    }[c]));
  }

  // Bold the (case-insensitive) occurrences of the typed text
  function highlight(text, q) {
    const safe = escapeHtml(text);
    if (!q) return safe;
    const re = new RegExp(escapeHtml(q).replace(/[.*+?^${}()|[\]\\]/g, '\\$&'), 'ig');
    return safe.replace(re, m => `<strong>${m}</strong>`);
  }

  function close() {
    list.hidden = true;
    list.innerHTML = '';
    items = [];
    active = -1;
    input.setAttribute('aria-expanded', 'false');
    input.removeAttribute('aria-activedescendant');
  }

  function setActive(i) {
    const options = list.children;
    if (active >= 0 && options[active]) options[active].classList.remove('active');
    active = i;
    if (active >= 0 && options[active]) {
      options[active].classList.add('active');
      options[active].scrollIntoView({ block: 'nearest' });
      input.setAttribute('aria-activedescendant', options[active].id);
    } else {
      input.removeAttribute('aria-activedescendant');
    }
  }

  function choose(i) {
    const item = items[i];
    if (!item) return;
    input.value = item.accession;
    close();
    input.form.requestSubmit ? input.form.requestSubmit() : input.form.submit();
  }

  function render(results, q) {
    items = results;
    active = -1;
    if (!results.length) { close(); return; }
    list.innerHTML = results.map((r, i) => `
      <li id="query-suggestion-${i}" role="option" data-index="${i}">
        <span class="sug-gene">${highlight(r.gene_name || '—', q)}</span>
        <span class="sug-name">${highlight(r.protein_name, q)}</span>
        <span class="sug-acc">${highlight(r.accession, q)}</span>
      </li>`).join('');
    list.hidden = false;
    input.setAttribute('aria-expanded', 'true');
  }

  function fetchSuggestions(q) {
    if (controller) controller.abort();
    controller = new AbortController();
    fetch(`${input.dataset.suggestUrl}?q=${encodeURIComponent(q)}`, { signal: controller.signal })
      .then(r => (r.ok ? r.json() : []))
      .then(results => {
        // Ignore a response for text the user has since changed
        if (input.value.trim() === q) render(results, q);
      })
      .catch(err => { if (err.name !== 'AbortError') close(); });
  }

  input.addEventListener('input', () => {
    clearTimeout(timer);
    const q = input.value.trim();
    if (q.length < MIN_CHARS) {
      if (controller) controller.abort();
      close();
      return;
    }
    timer = setTimeout(() => fetchSuggestions(q), DEBOUNCE_MS);
  });

  input.addEventListener('keydown', e => {
    if (list.hidden || !items.length) return;
    if (e.key === 'ArrowDown') {
      e.preventDefault();
      setActive((active + 1) % items.length);
    } else if (e.key === 'ArrowUp') {
      e.preventDefault();
      setActive(active <= 0 ? items.length - 1 : active - 1);
    } else if (e.key === 'Enter' && active >= 0) {
      e.preventDefault();
      choose(active);
    } else if (e.key === 'Escape') {
      close();
    }
  });

  // mousedown (not click) so the choice lands before the input's blur closes the list
  list.addEventListener('mousedown', e => {
    const li = e.target.closest('li[data-index]');
    if (!li) return;
    e.preventDefault();
    choose(Number(li.dataset.index));
  });
  list.addEventListener('mousemove', e => {
    const li = e.target.closest('li[data-index]');
    if (li && Number(li.dataset.index) !== active) setActive(Number(li.dataset.index));
  });

  input.addEventListener('blur', close);
})();
