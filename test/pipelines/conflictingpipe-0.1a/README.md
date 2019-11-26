It's possible to have multiple pipelines with the same name and version number. 
This shouldn't break any of the tools, but it will result in unpredictable 
behavior. Essentially, the first pipeline to match the name/version of the 
lookup will be returned. This selection should be consistent, but difficult to 
predict which will be found first, it comes down to the order of the search 
paths and the underlying filesystem.