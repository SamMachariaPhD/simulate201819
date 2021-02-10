#!/bin/bash
# Remove some files from remote repo but not on local repo
# Sam Macharia
git rm --cached *.csv # ignore all .csv files
git rm --cached *.txt
git rm --cached *.mp4
git rm --cached *.avi
git rm --cached *.png
git rm --cached *.html
git rm --cached *.vtk
git rm --cached *.svg
git rm --cached *.pdf
git rm --cached *.ipynb_checkpoints
