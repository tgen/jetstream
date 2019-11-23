# Tasks

Tasks are a set of key-value properties called directives. The task directives
determine _when_ and _how_ the task will be executed. But, some directives will 
be used differently depending on which backend you're running with.

Any additional directives can be added as well, they will not affect the runner
or execution of tasks. These will be stored in the workflow and could be useful 
for downstream meta-analysis of the projects themselves (eg. performance 
comparisons, or documenting the tasks)

Here is a list of the common task directives and how they're used: 

---

"Core" directives

- `name`

  A unique identifier for the task. If this is absent, it will be assigned 
  based on a hash of the task content. Naming tasks allows you to link 
  dependencies to this task by name. It's also used to determine the filename 
  for logs saved in projects.
  
- `cmd`

  Command to be executed by the backend, should be valid 
  [Bash](https://www.gnu.org/software/bash/).
  This is where the main "work" of the task should exist.

---

"Flow" directives 

Flow directives determine the order in which the tasks will be executed:

- `after` 

  This task will run after tasks named with each value. Supports sequences.
  
- `after-re`
  
  This task will run after tasks matching each given regex pattern. 
  Supports sequences.
  
- `before`

  This task will run before tasks named with each value. Supports sequences.

- `before-re`
  
  This task will run before tasks matching each given regex pattern. 
  Supports sequences.
  
- `input`:

  This task will run after tasks with output directives matching 
  each given value. Supports sequences.

- `input-re`:

  This task will run after tasks with output directives matching 
  each given regex pattern. Supports sequences.
  
- `output`: 

  This task will satisfy a matching `input` value requirement.
  Supports sequences.

---

"Execution" directives

- `exec`: 

  Python code that will execute in the runner process immediately 
  before sending the task to the backend. This feature can be used
  to modify the workflow (add tasks) while it's being run. Two local
  variables are added during execution: `task` and `runner`. The runner
  is important because it contains the current workflow. Any
  errors during execution will halt the runner immediately.
  
  Note: the workflow graph will always be recalculated after any exec 
  directive runs, and most work can be done with `cmd` directives. 

- `cpus`: 

  LocalBackend - Will reserve local cpus when launching cmd
  SlurmBackend - Passed as "-c" when requesting the job allocation
      
- `mem`: 

  SlurmBackend - Passed as "--mem" when requesting job allocation

- `reset`:

  When this task is reset, it will also reset other tasks in the workflow. 
  Typically when a task is reset, any descendants will also be reset. This 
  directive allows additional tasks to be reset. Options for this directive
  are: task name(s) or the special value "predecessors" which will reset any
  immediate upstream tasks. 
      
- `stdin`: 

  cmd stdin will be connected to this value

- `stdout`: 

  cmd stdout will be connected to this value

- `stderr`: 

  cmd stderr will be connected to this value


---

Any additional directives can be added as well, they will not affect the runner
or execution of tasks. These will be stored in the workflow and could be useful 
for downstream meta-analysis of the projects themselves (eg. performance 
comparisons, or documenting the tasks)

Here are some ideas that we've used before:

- `tags`: 
  
  Sequence of short descriptive tags that can be used for categorizing tasks
  for downstream analysis.
  
- `methods`: 

  Describe what this task does in plain language, later this can
  be used to generate a methods section for a project.
  
- `description`: 
  
  Describe this task for potential users
