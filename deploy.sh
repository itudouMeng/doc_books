tree -I ".git|mathjax">README.md
sed -i '1i # SUMMARY\
```' README.md
sed -i '$a ```' README.md

git add .
git commit -m "history"
git push origin master

docsify serve
