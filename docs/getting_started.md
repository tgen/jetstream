# Projects

Jetstream relies heavily on organizing work into projects. Most of the commands covered in this guide are going to require that you are working inside of a project directory. To get started

# Pipelines

Jetstream runs computational pipelines by reading workflow templates. Some templates are included with the package and can be run 

Pipeline are a set of tasks needed to complete a computational goal. They could be as simple as a list of sequential commands to run, or a complicated network of modeling dependency and concurrency. 


# Loading templates


## Altering the search path preference

```
[~/myproject]$ jetstream_pipelines hello_world.yaml
Loaded template: templates/hello_world.yaml
...
```

```
[rrichholt@it5687:~/myproject]$ jetstream_pipelines --no-project-templates hello_world.yaml
Loaded template: /Users/rrichholt/jetstream/jetstream_pipelines/templates/hello_world.yaml
...
```
