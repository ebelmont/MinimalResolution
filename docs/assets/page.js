(function () {
  "use strict";

  // Resolve this script's own directory so mermaid.min.js (shipped alongside
  // it, no CDN) can be loaded with the right relative path regardless of how
  // deep the current page is nested (docs/ vs docs/pipelines/).
  var thisScript = document.currentScript;
  var assetsBase = thisScript ? thisScript.src.replace(/page\.js(\?.*)?$/, "") : "";

  // Theme toggle: explicit choice overrides prefers-color-scheme, persisted per-browser.
  var root = document.documentElement;
  var stored = null;
  try { stored = localStorage.getItem("mr-docs-theme"); } catch (e) {}
  if (stored === "light" || stored === "dark") root.setAttribute("data-theme", stored);

  function currentTheme() {
    var attr = root.getAttribute("data-theme");
    if (attr) return attr;
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }

  function initThemeToggle() {
    var btn = document.querySelector("[data-theme-toggle]");
    if (!btn) return;
    var setLabel = function () { btn.textContent = currentTheme() === "dark" ? "☀ light" : "● dark"; };
    setLabel();
    btn.addEventListener("click", function () {
      var next = currentTheme() === "dark" ? "light" : "dark";
      root.setAttribute("data-theme", next);
      try { localStorage.setItem("mr-docs-theme", next); } catch (e) {}
      setLabel();
    });
  }

  // Build a right-rail table of contents from h2/h3 headings, only if there are enough to matter.
  function initToc() {
    var toc = document.getElementById("toc");
    var article = document.querySelector("article");
    if (!toc || !article) return;
    var headings = article.querySelectorAll("h2[id], h3[id]");
    if (headings.length < 4) { toc.remove(); return; }
    var list = document.createElement("ul");
    headings.forEach(function (h) {
      var li = document.createElement("li");
      if (h.tagName === "H3") li.className = "lvl-3";
      var a = document.createElement("a");
      a.href = "#" + h.id;
      a.textContent = h.textContent.replace(/\s*#$/, "");
      li.appendChild(a);
      list.appendChild(li);
    });
    var title = document.createElement("div");
    title.className = "toc-title";
    title.textContent = "On this page";
    toc.appendChild(title);
    toc.appendChild(list);

    var links = list.querySelectorAll("a");
    var byId = {};
    links.forEach(function (a) { byId[a.getAttribute("href").slice(1)] = a; });
    var obs = new IntersectionObserver(function (entries) {
      entries.forEach(function (entry) {
        var link = byId[entry.target.id];
        if (!link) return;
        link.classList.toggle("active", entry.isIntersecting);
      });
    }, { rootMargin: "0px 0px -70% 0px" });
    headings.forEach(function (h) { obs.observe(h); });
  }

  // Mermaid diagrams: progressive enhancement, loaded from a local copy
  // shipped in docs/assets/ (no CDN -- works offline and works when this
  // page is hosted standalone). If the script fails to load for any reason,
  // the preformatted source text in .mermaid-wrap is already legible on its own.
  function initMermaid() {
    var nodes = document.querySelectorAll(".mermaid");
    if (!nodes.length) return;
    var script = document.createElement("script");
    script.src = assetsBase + "mermaid.min.js";
    script.onload = function () {
      var dark = document.documentElement.getAttribute("data-theme") === "dark" ||
          (!document.documentElement.getAttribute("data-theme") && window.matchMedia("(prefers-color-scheme: dark)").matches);
      window.mermaid.initialize({ startOnLoad: false, theme: dark ? "dark" : "default" });
      window.mermaid.run({ querySelector: ".mermaid" }).catch(function () {});
    };
    document.body.appendChild(script);
  }

  document.addEventListener("DOMContentLoaded", function () {
    initThemeToggle();
    initToc();
    initMermaid();
  });
})();
