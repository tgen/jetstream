# TL;DR

- Jetstream includes:
	- Comprehensive genomics pipelines (Specifically tailored for the TGen infrastructure)
	- Core workflow engine that runs these pipelines and allows users to build new pipelines
	- A handful of additional scripts and Python modules for various genomics tasks

- `jetstream_pipelines` is used to run pipelines (included and custom)

- Jetstream pipelines have to be run on a project directory, and these have some simple, but essential requirements.

- `jetstream` is used for creating project directories

- Jetstream pipelines are stored as yaml-jinja text documents called [templates](templates).


# Intro

# note: often times when I read documents for a project, I find myself thinking that they dont spend enough time introducing the project at the basic/beginner level. I want this section to be a very clear introduction to the application, but ive been working so closely on it that its difficult for me to step back and see it from an outside perspective.


# Projects

Jetstream relies heavily on organizing work as projects. Most of the commands covered in this guide are going to require that you are working inside of a project, and all the built-in pipelines assume that you're inside of a project. This section covers a description of the Jetstream project and the rationale behind the structure. 

Jetstream Projects are just directories on a filesystem that include a few simple, but essential requirements. You can 


# Pipelines

Jetstream runs computational pipelines by combining project data with workflow templates. Some templates are included with the package and can be run 

Pipeline are a set of tasks needed to complete a computational goal. They could be as simple as a list of sequential commands to run, or a complicated network of modeling dependency and concurrency. 


# Templates

> This section describes the underlying workflow engine and how to create new pipelines. If you just want to use the built-in pipelines, you can skip this part for now.

The Jetstream workflow engine establishes some ground rules for describing commands, and provides a set of tools for accessing project data. Beyond this, the pipelines are completely user customizeable and separate from the application code. This allows rapid development of the built-in pipelines, and also creates a canvas for users to architect pipelines that suites their needs. 

Templates are awesome. They allow us to generate very complex workflows with a minimalistic syntax. They also make workflows more transparent - the flow control, or order in which the tasks will be completed, is _entirely_ expressed by the template. There are no internal, engine-level, workflow procedures. After the template is rendered, it is an exact description of the commands that will be executed. 

(# TODO make graphic)
```

Code/Scripts 		<----------|Flexibility|----------->     CLI/GUI Applications
 
 								 Templates

+ Flexibilty              + Very flexible                    + No user development time needed 
- Speed					  + Limited/Predictable behavior     + Very predicatble results
- Predictability          - Small learning curve             - Flexibile only by adding more options

 ```


https://en.wikipedia.org/wiki/Template_processor


http://jinja.pocoo.org/docs/2.10/templates/


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
