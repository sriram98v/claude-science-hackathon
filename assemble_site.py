"""Assemble the Zola site inputs from the portable blog contract (blog/ only).

Reads ONLY blog/content.md, blog/media/, and blog/full.html -- never the notebook --
and materialises them into site/ so `zola build` can run:

  blog/media/            -> site/static/media/
  blog/full.html         -> site/static/full.html
  blog/content.md        -> site/content/_index.md   (two-column reader, template=reader.html)
                         -> site/content/article.md  (single-column,   template=article.html)

This is the single place the two Zola pages get their (identical) body, so there is
still one source of truth. Both CI (deploy-web) and local verification call this, then
run `zola build` in site/. Run:  python assemble_site.py
"""
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
BLOG = os.path.join(HERE, "blog")
SITE = os.path.join(HERE, "site")

READER_FM = "+++\ntemplate = \"reader.html\"\n+++\n\n"
ARTICLE_FM = "+++\ntitle = \"Article\"\npath = \"article\"\ntemplate = \"article.html\"\n+++\n\n"


def main():
    content_md = os.path.join(BLOG, "content.md")
    media_src = os.path.join(BLOG, "media")
    full_src = os.path.join(BLOG, "full.html")
    for p in (content_md, media_src):
        if not os.path.exists(p):
            raise SystemExit(f"missing blog input: {p} -- run export_blog.py first")

    with open(content_md, encoding="utf-8") as f:
        body = f.read()

    # static/media (mirror, replacing any stale copy)
    media_dst = os.path.join(SITE, "static", "media")
    if os.path.isdir(media_dst):
        shutil.rmtree(media_dst)
    shutil.copytree(media_src, media_dst)

    # full.html reference page
    os.makedirs(os.path.join(SITE, "static"), exist_ok=True)
    if os.path.exists(full_src):
        shutil.copy2(full_src, os.path.join(SITE, "static", "full.html"))

    # two content files, same body, different front matter
    content_dir = os.path.join(SITE, "content")
    os.makedirs(content_dir, exist_ok=True)
    with open(os.path.join(content_dir, "_index.md"), "w", encoding="utf-8") as f:
        f.write(READER_FM + body)
    with open(os.path.join(content_dir, "article.md"), "w", encoding="utf-8") as f:
        f.write(ARTICLE_FM + body)

    n_media = sum(len(files) for _, _, files in os.walk(media_dst))
    print(f"[assemble_site] wrote site/content/_index.md + article.md, "
          f"site/static/media ({n_media} files), site/static/full.html")


if __name__ == "__main__":
    main()
