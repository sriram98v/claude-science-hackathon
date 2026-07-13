(function () {
  var root = document.documentElement;
  var wide = window.matchMedia("(min-width: 1400px)");
  var textCol = document.getElementById("col-text");
  var mediaCol = document.getElementById("col-media");
  var reader = document.querySelector(".reader");
  var textBody = document.querySelector(".text-body");

  var anchors = Array.prototype.slice.call(document.querySelectorAll(".media-anchor"));
  var blocks = {};
  anchors.forEach(function (a) {
    blocks[a.dataset.media] = document.getElementById("media-" + a.dataset.media);
  });

  /* ---- move media blocks between the right column and inline ---- */
  function toColumn() {
    anchors.forEach(function (a) {
      var b = blocks[a.dataset.media];
      if (b && b.parentNode !== mediaCol) mediaCol.appendChild(b);
    });
  }
  function toInline() {
    anchors.forEach(function (a) {
      var b = blocks[a.dataset.media];
      if (b && a.nextSibling !== b) a.parentNode.insertBefore(b, a.nextSibling);
    });
  }

  /* tag blocks whose table is wider than a portrait A4 text column (~640px) so
     print CSS can send them to their own landscape page */
  function tagWideTables() {
    anchors.forEach(function (a) {
      var b = blocks[a.dataset.media];
      if (!b) return;
      var wideT = false;
      Array.prototype.forEach.call(b.querySelectorAll("table"), function (t) {
        if (Math.max(t.scrollWidth, t.getBoundingClientRect().width) > 640) wideT = true;
      });
      b.classList.toggle("wide-table", wideT);
    });
  }

  /* ---- current scroll container depends on layout ---- */
  function scroller() { return wide.matches ? textCol : reader; }

  /* ---- build the index panel ---- */
  var list = document.getElementById("toc-list");
  var heads = textBody ? Array.prototype.slice.call(textBody.querySelectorAll("h2, h3")) : [];
  var tocById = {};
  heads.forEach(function (h) {
    if (!h.id) h.id = h.textContent.trim().replace(/\s+/g, "-");
    var a = document.createElement("a");
    a.href = "#" + h.id;
    a.textContent = h.textContent.replace(/¶/g, "").trim();
    a.className = h.tagName.toLowerCase() === "h3" ? "toc-sub" : "toc-top";
    a.dataset.target = h.id;
    list.appendChild(a);
    tocById[h.id] = a;
    a.addEventListener("click", function (e) {
      e.preventDefault();
      jumpTo(h);
      if (!wide.matches) document.body.classList.remove("toc-open");
    });
  });

  /* ---- media sync (left drives right; suppressed while hovering media) ---- */
  var overMedia = false;
  mediaCol.classList.add("focus-mode");   // dim non-active blocks by default
  mediaCol.addEventListener("mouseenter", function () {
    overMedia = true; mediaCol.classList.remove("focus-mode");   // full opacity while inspecting
  });
  mediaCol.addEventListener("mouseleave", function () {
    overMedia = false; mediaCol.classList.add("focus-mode");
  });

  function activateMedia(k, doScroll) {
    Object.keys(blocks).forEach(function (key) {
      if (blocks[key]) blocks[key].classList.toggle("active", key === k);
    });
    if (doScroll && wide.matches && !overMedia && blocks[k]) {
      var blk = blocks[k];
      var target = blk.offsetTop - (mediaCol.clientHeight - blk.offsetHeight) / 2;
      var max = mediaCol.scrollHeight - mediaCol.clientHeight;
      mediaCol.scrollTo({ top: Math.max(0, Math.min(target, max)), behavior: "smooth" });
    }
  }

  var syncObs = null, spyObs = null;
  function setupObservers() {
    if (syncObs) syncObs.disconnect();
    if (spyObs) spyObs.disconnect();
    var rootEl = scroller();

    /* right column follows the anchor near the top of the reading area */
    if (wide.matches) {
      syncObs = new IntersectionObserver(function (entries) {
        entries.forEach(function (e) {
          if (e.isIntersecting) activateMedia(e.target.dataset.media, true);
        });
      }, { root: rootEl, rootMargin: "-30% 0px -65% 0px", threshold: 0 });
      anchors.forEach(function (a) { syncObs.observe(a); });
    }

    /* index highlight */
    spyObs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) {
          heads.forEach(function (h) { if (tocById[h.id]) tocById[h.id].classList.remove("active"); });
          if (tocById[e.target.id]) tocById[e.target.id].classList.add("active");
        }
      });
    }, { root: rootEl, rootMargin: "-12% 0px -78% 0px", threshold: 0 });
    heads.forEach(function (h) { spyObs.observe(h); });
  }

  /* ---- jump to a heading and sync the right column ---- */
  function jumpTo(h) {
    h.scrollIntoView({ behavior: "smooth", block: "start" });
    /* find the first media anchor at or after this heading */
    var k = null;
    for (var i = 0; i < anchors.length; i++) {
      var rel = h.compareDocumentPosition(anchors[i]);
      if (rel & Node.DOCUMENT_POSITION_FOLLOWING) { k = anchors[i].dataset.media; break; }
    }
    if (k) activateMedia(k, true);
  }

  /* ---- apply layout for the current breakpoint ---- */
  function applyLayout() {
    if (wide.matches) { toColumn(); } else { toInline(); }
    setupObservers();
    tagWideTables();
  }
  applyLayout();
  if (wide.addEventListener) wide.addEventListener("change", applyLayout);
  else wide.addListener(applyLayout);

  /* hamburger toggles ONLY the index: collapse it (wide) or slide overlay (narrow) */
  document.getElementById("toc-toggle").addEventListener("click", function () {
    if (wide.matches) document.body.classList.toggle("toc-collapsed");
    else document.body.classList.toggle("toc-open");
  });

  /* ---- dark-mode toggle ---- */
  var themeBtn = document.getElementById("theme-toggle");
  function syncIcon() {
    themeBtn.innerHTML = root.getAttribute("data-theme") === "dark" ? "&#9728;" : "&#9790;";
  }
  themeBtn.addEventListener("click", function () {
    var next = root.getAttribute("data-theme") === "dark" ? "light" : "dark";
    root.setAttribute("data-theme", next);
    try { localStorage.setItem("blog-theme", next); } catch (e) {}
    syncIcon();
  });
  syncIcon();

  /* ---- scroll-to-top (acts on the active scroller) ---- */
  var stt = document.getElementById("scroll-top");
  function onScroll() { stt.classList.toggle("show", scroller().scrollTop > 400); }
  textCol.addEventListener("scroll", onScroll);
  reader.addEventListener("scroll", onScroll);
  stt.addEventListener("click", function () {
    scroller().scrollTo({ top: 0, behavior: "smooth" });
    if (wide.matches) mediaCol.scrollTo({ top: 0, behavior: "smooth" });
  });

  /* ---- export / print: single-column, media inline; restore afterward ---- */
  document.getElementById("export-pdf").addEventListener("click", function () { window.print(); });
  window.addEventListener("beforeprint", function () { toInline(); tagWideTables(); });
  window.addEventListener("afterprint", applyLayout);

  /* ---- click a figure to zoom (single-column mode only) ---- */
  var lb = document.getElementById("lightbox");
  var lbImg = lb.querySelector("img");
  function closeLightbox() { lb.classList.remove("open"); }
  lb.addEventListener("click", closeLightbox);
  document.addEventListener("keydown", function (e) { if (e.key === "Escape") closeLightbox(); });
  Array.prototype.forEach.call(document.querySelectorAll(".media-block img"), function (img) {
    img.addEventListener("click", function (e) {
      if (wide.matches) return;          // only when stacked in a single column
      e.stopPropagation();
      lbImg.src = img.currentSrc || img.src;
      lb.classList.add("open");
    });
  });
})();
