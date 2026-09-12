/* Builds the table of contents for the current page from its headings, keeps
   the entry for the section being read highlighted, and drives the sidebar
   drawer on narrow screens. */
(function () {
  'use strict';

  var header = document.querySelector('.site-header');
  var sidebar = document.getElementById('sidebar');
  var scrim = document.querySelector('.nav-scrim');

  /* ------------------------------------------------------ drawer toggle */

  var toggle = document.querySelector('.nav-toggle');

  function setNav(open) {
    document.body.classList.toggle('nav-open', open);
    if (toggle) toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
  }

  if (toggle) {
    toggle.addEventListener('click', function () {
      setNav(!document.body.classList.contains('nav-open'));
    });
  }
  if (scrim) scrim.addEventListener('click', function () { setNav(false); });
  document.addEventListener('keydown', function (e) {
    if (e.key === 'Escape') setNav(false);
  });

  /* -------------------------------------------------- build the contents */

  var main = document.getElementById('content');
  var slot = document.getElementById('page-toc');
  if (!main || !slot || !sidebar) return;

  /* Heading text, with any mathematics reduced to its TeX source.  Several
     headings are partly or wholly mathematics, so dropping it outright would
     leave entries like "and , as discretised".

     The maths can be in either of two states, because this script and MathJax
     both load asynchronously.  Untouched, it is plain text between \( \)
     delimiters, exactly as kramdown wrote it.  Once MathJax has run, each
     formula is present three times over -- a hidden preview, the rendered
     output, and the original TeX in a <script> -- and the class of the
     rendered container depends on the output processor (MathJax_CHTML here),
     so match on the prefix rather than on one class. */
  function headingText(h) {
    var clone = h.cloneNode(true);

    var classed = clone.querySelectorAll('[class]');
    for (var i = classed.length - 1; i >= 0; i--) {
      var cls = classed[i].getAttribute('class') || '';
      if (cls.indexOf('MathJax') > -1 && classed[i].parentNode) {
        classed[i].parentNode.removeChild(classed[i]);
      }
    }

    var math = clone.querySelectorAll('script[type^="math/tex"]');
    for (var j = 0; j < math.length; j++) {
      math[j].parentNode.replaceChild(document.createTextNode(math[j].textContent), math[j]);
    }

    return clone.textContent
      .replace(/\\[()[\]]/g, ' ')   /* \( \) \[ \] */
      .replace(/\$\$?/g, ' ')       /* $ and $$ */
      .replace(/\\[,!;:>]/g, '')    /* thin spaces and friends */
      .replace(/[{}\\]/g, '')       /* braces, and the slash of \alpha etc. */
      .replace(/\s+/g, ' ')
      .replace(/\s+([,.;:)])/g, '$1')  /* the delimiters left gaps behind */
      .trim();
  }

  var headings = main.querySelectorAll('h2, h3');
  if (headings.length < 2) return;

  var list = document.createElement('ul');
  list.className = 'toc';
  var links = [];
  var targets = [];

  for (var k = 0; k < headings.length; k++) {
    var h = headings[k];
    var text = headingText(h);
    if (!text) continue;
    if (!h.id) h.id = 'section-' + (k + 1);

    var li = document.createElement('li');
    li.className = 'toc-' + h.tagName.toLowerCase();
    var a = document.createElement('a');
    a.href = '#' + h.id;
    a.textContent = text;
    li.appendChild(a);
    list.appendChild(li);
    links.push(a);
    targets.push(h);
  }

  if (!links.length) return;
  slot.appendChild(list);

  /* Tapping an entry on a narrow screen should reveal the section, not leave
     the drawer covering it. */
  list.addEventListener('click', function (e) {
    if (e.target.tagName === 'A') setNav(false);
  });

  /* ------------------------------------------------------- reading position */

  var current = -1;

  function scrollIntoSidebar(a) {
    var r = a.getBoundingClientRect();
    var s = sidebar.getBoundingClientRect();
    if (r.top < s.top + 8) sidebar.scrollTop += r.top - s.top - 8;
    else if (r.bottom > s.bottom - 8) sidebar.scrollTop += r.bottom - s.bottom + 8;
  }

  function update() {
    var cut = (header ? header.offsetHeight : 0) + 24;
    /* -1 while still in the material above the first heading. */
    var index = -1;
    for (var i = 0; i < targets.length; i++) {
      if (targets[i].getBoundingClientRect().top - cut <= 0) index = i;
      else break;
    }
    /* At the very bottom the last section may never clear the cut-off. */
    if (window.innerHeight + window.pageYOffset >= document.body.offsetHeight - 2) {
      index = targets.length - 1;
    }
    if (index === current) return;
    if (current >= 0) links[current].classList.remove('is-active');
    current = index;
    if (current >= 0) {
      links[current].classList.add('is-active');
      scrollIntoSidebar(links[current]);
    }
  }

  var queued = false;
  function onScroll() {
    if (queued) return;
    queued = true;
    window.requestAnimationFrame(function () {
      queued = false;
      update();
    });
  }

  window.addEventListener('scroll', onScroll, { passive: true });
  window.addEventListener('resize', onScroll);
  window.addEventListener('load', function () {
    onScroll();
    /* Typesetting moves everything down, so measure again once it is done. */
    if (window.MathJax && window.MathJax.Hub && window.MathJax.Hub.Queue) {
      window.MathJax.Hub.Queue(onScroll);
    }
  });
  update();
})();
