Extremly short introduction to shell
====================================

.. highlight:: bash

.. attention::

    If you are already (to any extent) familiar with Linux shell,
    skip directly to :ref:`Task 2 <bash-task2>`.

Shell is the direct interface between the user and the system kernel.

.. note::

    Bash is one of several shell interpreters. Other shells you can meet
    are: sh, csh, tcsh, zsh, ash, etc...  The bash shell is the default
    standad shell on most GNU/linux based systems and basic commands are
    however almost same across all the shell variants.

.. admonition:: Task 1

    Open any terminal application - gnome-terminal, kterm, xterm.
    On many Lunux distributions you can achieve that by pushing ``CTRL+ALT+T``.
    You will, most probably, be greeted by Bash prompt::

       user@machine:~$


Basic shell commands
--------------------

``whoami``
      who am I?
``pwd``
      print working directory
``ls``
      list working directory
``ls -l``
      detailed list of working directory
``cd``
      change directory
``cp``
      copy files
``mv``
      move files
``rm``
      delete file
``mkdir``
      create directory
``rm -rf``
      remove directory recursively (including its contents)

If the command is not recognized as shell command, then the folders in
the enviroment variable ``PATH`` are searched for an executable with
that name.


.. hint::

    Help can be obtained using commands ``man``, ``info``,
    ``apropos``. Try for instance::

        man ls
        apropos math

    To exit a ``man`` page hit ``q`` key.


Remote connection - ssh
-----------------------

.. _bash-task2:

.. attention::

    If you are already familiar with SSH login to TUC machines,
    skip to :ref:`Python tutorial <python-intro>`.

.. admonition:: Task 2

    Use ``ssh`` utility to connect to a remote system::

        ssh -X -C tyche

    .. note::

        * ``-X`` enable ``X11`` forwarding (allows processes on the remote machine
          opening windows of graphical applications on the local machine)
        * ``-C`` enables compression which is mainly beneficial for access from
          a remote network
        * ``tyche`` stands here for machine ``tyche.mathematik.tu-chemnitz.de``;
          username on the local machine is used by default to login to the remote
          machine; the machine is not accessible from outside the univerity, so
          one would login through a jump host::

            helmut@local_machine:~$ ssh -X -C user@login.tu-chemnitz.de
            user@login:~$ ssh -X -C tyche
            user@tyche:~$

          Alternatively one can use VPN.


Sending signals and process management
--------------------------------------

Every running program can be in the state *running* or *suspended*.
With respect to a shell it can be in the *foreground* or *background*.

``jobs``
     list jobs started from given shell with their shell ids
``fg``
     send job to foreground
``bg``
     send job to background, this is equivalent to starting the job by ``job_command &``

Every well behaved job listens to signals.

``ps``
    list all your jobs with their process ids
``kill jobid``
    send *terminate* signal to a job, job id can be its process id or its shell id (given by ``%<num>``)
``kill -9 jobid``
    send *kill* signal (immediate request to end) to a job

Pressing following control keys will send signal to the foreground job

``^C``
           terminate signal
``^Z``
           suspend signal
``^D``
           end of input signal
