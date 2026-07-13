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
  var heads = textBody ? Array.prototype.slice.call(textBody.querySelectorAll("h1, h2")) : [];
  /* the first h1 is the article title itself, not an index entry */
  if (heads.length && heads[0].tagName.toLowerCase() === "h1") heads.shift();
  var tocById = {};
  heads.forEach(function (h) {
    h.id = h.textContent.replace(/¶/g, "").trim().replace(/\s+/g, "-");
    var a = document.createElement("a");
    a.href = "#" + h.id;
    a.textContent = h.textContent.replace(/¶/g, "").trim();
    a.className = h.tagName.toLowerCase() === "h2" ? "toc-sub" : "toc-top";
    list.appendChild(a);
    tocById[h.id] = a;
    a.addEventListener("click", function () {
      /* only the overlay needs dismissing; the wide panel stays put */
      if (!wide.matches) document.body.classList.remove("toc-open");
    });
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
