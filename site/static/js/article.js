/* Single-column article: one text column with the Contents panel pinned on the left on
   wide screens (like the two-column reader) and a slide-in overlay on narrow screens.
   No media column -- figures/tables stay inline at their anchors and zoom on click. */
(function () {
  var root = document.documentElement;
  var wide = window.matchMedia("(min-width: 1400px)");
  var textCol = document.getElementById("col-text");
  var reader = document.querySelector(".reader");
  var textBody = document.querySelector(".text-body");

  /* the scroll container depends on layout: on wide the text column scrolls inside the
     grid; on narrow the whole reader scrolls as one. */
  function scroller() { return wide.matches ? textCol : reader; }

  /* ---- build the index panel ---- */
  var list = document.getElementById("toc-list");
  var heads = textBody ? Array.prototype.slice.call(textBody.querySelectorAll("h1, h2, h3, h4")) : [];
  /* the first h1 is the article title itself, not an index entry */
  if (heads.length && heads[0].tagName.toLowerCase() === "h1") heads.shift();
  var tocById = {};
  /* Build a collapsible tree: h2 sections hold h3 subsections hold h4 items. A caret on
     any node with children toggles its subtree; the label itself still navigates (native
     anchor jump). */
  var tocStack = [{ level: 1, container: list }];
  heads.forEach(function (h) {
    h.id = h.textContent.replace(/¶/g, "").trim().replace(/\s+/g, "-");
    var level = parseInt(h.tagName.substring(1), 10);   // 2, 3, 4
    while (tocStack.length > 1 && tocStack[tocStack.length - 1].level >= level) tocStack.pop();
    var parentBox = tocStack[tocStack.length - 1].container;

    var item = document.createElement("div");
    item.className = "toc-item toc-lvl" + level;
    var row = document.createElement("div");
    row.className = "toc-row";
    var toggle = document.createElement("button");
    toggle.className = "toc-toggle";
    toggle.setAttribute("aria-label", "Collapse");
    toggle.setAttribute("aria-expanded", "true");
    var a = document.createElement("a");
    a.href = "#" + h.id;
    a.textContent = h.textContent.replace(/¶/g, "").trim();
    a.className = "toc-link toc-h" + level;
    row.appendChild(toggle);
    row.appendChild(a);
    item.appendChild(row);
    var kids = document.createElement("div");
    kids.className = "toc-children";
    item.appendChild(kids);
    parentBox.appendChild(item);

    tocById[h.id] = a;
    a.addEventListener("click", function () {
      /* only the overlay needs dismissing; the wide panel stays put */
      if (!wide.matches) document.body.classList.remove("toc-open");
    });
    toggle.addEventListener("click", function (e) {
      e.preventDefault();
      e.stopPropagation();
      var collapsed = item.classList.toggle("collapsed");
      toggle.setAttribute("aria-expanded", collapsed ? "false" : "true");
      toggle.setAttribute("aria-label", collapsed ? "Expand" : "Collapse");
    });

    tocStack.push({ level: level, container: kids });
  });
  /* leaves drop their caret (kept as blank space); parents start collapsed */
  Array.prototype.forEach.call(list.querySelectorAll(".toc-item"), function (it) {
    var box = it.querySelector(":scope > .toc-children");
    if (!box || box.children.length === 0) {
      it.classList.add("toc-leaf");
    } else {
      it.classList.add("collapsed");
      var tg = it.querySelector(":scope > .toc-row > .toc-toggle");
      if (tg) { tg.setAttribute("aria-expanded", "false"); tg.setAttribute("aria-label", "Expand"); }
    }
  });

  /* ---- index highlight, rooted on the active scroll container ---- */
  var spyObs = null;
  function setupSpy() {
    if (spyObs) spyObs.disconnect();
    spyObs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) {
          heads.forEach(function (h) { if (tocById[h.id]) tocById[h.id].classList.remove("active"); });
          if (tocById[e.target.id]) tocById[e.target.id].classList.add("active");
        }
      });
    }, { root: scroller(), rootMargin: "-12% 0px -78% 0px", threshold: 0 });
    heads.forEach(function (h) { spyObs.observe(h); });
  }
  setupSpy();
  if (wide.addEventListener) wide.addEventListener("change", setupSpy);
  else wide.addListener(setupSpy);

  /* hamburger: collapse the pinned panel (wide) or slide the overlay (narrow) */
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
  });

  /* ---- export / print ---- */
  document.getElementById("export-pdf").addEventListener("click", function () { window.print(); });

  /* ---- click a figure to zoom (media is always inline here, so at all widths) ---- */
  var lb = document.getElementById("lightbox");
  var lbImg = lb.querySelector("img");
  function closeLightbox() { lb.classList.remove("open"); }
  lb.addEventListener("click", closeLightbox);
  document.addEventListener("keydown", function (e) { if (e.key === "Escape") closeLightbox(); });
  Array.prototype.forEach.call(document.querySelectorAll(".media-block img"), function (img) {
    img.addEventListener("click", function (e) {
      e.stopPropagation();
      lbImg.src = img.currentSrc || img.src;
      lb.classList.add("open");
    });
  });
})();
