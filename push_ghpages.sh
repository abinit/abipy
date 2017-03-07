#!/usr/bin/env bash
# http://www.willmcginnis.com/2016/02/29/automating-documentation-workflow-with-sphinx-and-github-pages/

HTML_PATH=./docs/_build/html

# build the docs
cd docs
make clean
make html
cd ..
# commit and push
git add -A
git commit -m "building and pushing docs"
git push origin master
# switch branches and pull the data we want
git checkout gh-pages
rm -rf .
touch .nojekyll
#git checkout master docs/build/html
#mv ./docs/build/html/* ./
git checkout master ${SITE_PATH}
mv ${SITE_PATH}/html/* ./
rm -rf ./docs
git add -A
git commit -m "publishing updated docs..."
git push origin gh-pages
# switch back
git checkout master