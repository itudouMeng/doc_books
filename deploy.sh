tree -I ".git|mathjax">README.md
sed -i '1i # SUMMARY\
```' README.md
sed -i '$a ```' README.md

git add .
<<<<<<< HEAD
git commit -m "maozedong"
git push
=======
git commit -m "history"
git push origin master
>>>>>>> 0cecec68ab8d7d5bc5e90da80e1f762a74154481

docsify serve
