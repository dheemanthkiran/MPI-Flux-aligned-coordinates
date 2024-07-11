# COLLECTION OF PYTHON TOOLS FOR ANALYTICAL MHD EQUILIBRIA 

## ipython notebooks

Please before commiting ipython notebooks `*.ipynb`, strip the output and metadata from file, by configuring once inside the local repo (requirement is jupyter nbconvert):
```
git config --local filter.strip-notebook-output.clean 'jupyter-nbconvert  --ClearOutputPreprocessor.enabled=True --ClearMetadataPreprocessor.enabled=True  --to=notebook --stdin --stdout --log-level=ERROR'
```

