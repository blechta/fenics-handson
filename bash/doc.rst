Extremly short introduction to bash
===================================

Bash is one of several shell interpreters. Other shells you can meet
are: sh, csh, tcsh, zsh, ash, etc...  The bash shell is the default
standad shell on most gnu/linux based systems and basic commands are
however almost same across all the shell variants. 

Shell is the direct interface between the user and the system kernel.


   **Task 1.** Open any terminal application - gnome-terminal, kterm, xterm.
   You will, most probably, be greeted by bash prompt

   .. code-block:: bash

      r0:~> 

Basic shell commands
--------------------

``whoami`` 
      kdo jsem?

``pwd``
      print working directory

``ls``
      list working directory

``ls -l``
      detailed list of working directory

``cd``
      change directory 

``cp``


``mv``



If the command is not recognized as shell command, then the folders in
the enviroment variable ``PATH`` are searched for an executable with
that name.


Remote connection - ssh
-----------------------

some text

   **Task 2.** Use ``ssh`` utility to connect to a remote system.
 
   .. code-block: bash

   machine:~> ssh -C -Y your_login@r0.karlin.mff.cuni.cz



Sending signals and process management
--------------------------------------

Every running program can be in the state *running* or *suspended*.
With respect to a shell it can be in the *foreground* or *background*.

| ``jobs``    list jobs started from given shell with their shell ids
| ``fg``      send job to foreground
| ``bg``      send job to background, this is equivalent to starting the job by ``job_command &``

Every well behaved job listens to signals.

| ``ps``                 list all your jobs with their process ids
| ``kill jobid``         send terminate signal to a job, job id can be its process id or its shell id (including %)
| ``kill -9 jobid``      change directory 

Pressing following control keys will send signal to the foreground job

| ``^C``           terminate signal
| ``^Z``           suspend signal
| ``^D``           end of input signal

Cluster job queues usage
----------------------------------

.. todo::

   Introduction to classes is needed to understand DOLFIN code.
