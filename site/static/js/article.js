/* Single-column article: no media column, no left/right sync. Media stays inline at
   its anchor; the TOC is a slide-in overlay at every width; figures zoom on click. */
(function () {
  var root = document.documentElement;
  var reader = document.querySelector(".reader");
  var textBody = document.querySelector(".text-body");

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
    list.appendChild(a);
    tocById[h.id] = a;
    a.addEventListener("click", function () {
      document.body.classList.remove("toc-open");
    });
  });

  /* index highlight */
  var spyObs = new IntersectionObserver(function (entries) {
    entries.forEach(function (e) {
      if (e.isIntersecting) {
        heads.forEach(function (h) { if (tocById[h.id]) tocById[h.id].classList.remove("active"); });
        if (tocById[e.target.id]) tocById[e.target.id].classList.add("active");
      }
    });
  }, { root: reader, rootMargin: "-12% 0px -78% 0px", threshold: 0 });
  heads.forEach(function (h) { spyObs.observe(h); });

  /* hamburger toggles the index overlay */
  document.getElementById("toc-toggle").addEventListener("click", function () {
    document.body.classList.toggle("toc-open");
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

  /* ---- scroll-to-top ---- */
  var stt = document.getElementById("scroll-top");
  reader.addEventListener("scroll", function () {
    stt.classList.toggle("show", reader.scrollTop > 400);
  });
  stt.addEventListener("click", function () {
    reader.scrollTo({ top: 0, behavior: "smooth" });
  });

  /* ---- export / print ---- */
  document.getElementById("export-pdf").addEventListener("click", function () { window.print(); });

  /* ---- click a figure to zoom (all widths -- always single column) ---- */
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
