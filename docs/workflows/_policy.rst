.. _task_policy:

^^^^^^^^^^
TaskPolicy
^^^^^^^^^^

At this point, you may wonder why we need to specify all these parameters in the configuration file.
The reason is that, before submitting a job to a resource manager, AbiPy will use the autoparal 
feature of ABINIT to get all the possible parallel configurations with ``ncpus <= max_cores``. 
On the basis of these results, AbiPy selects the "optimal" one, and changes the ABINIT input file 
and the submission script accordingly .
(this is a very useful feature, especially for calculations done with ``paral_kgb=1`` that require 
the specification of ``npkpt``, ``npfft``, ``npband``, etc).
If more than one `QueueAdapter` is specified, AbiPy will first compute all the possible 
configuration and then select the "optimal" `QueueAdapter` according to some kind of policy

In some cases, you may want to enforce some constraint on the "optimal" configuration. 
For example, you may want to select only those configurations whose parallel efficiency is greater than 0.7 
and whose number of MPI nodes is divisible by 4. 
One can easily enforce this constraint via the ``condition`` dictionary whose syntax is similar to 
the one used in ``mongodb``

.. code-block:: yaml

    policy:
        autoparal: 1
        max_ncpus: 10
        condition: {$and: [ {"efficiency": {$gt: 0.7}}, {"tot_ncpus": {$divisible: 4}} ]}

The parallel efficiency is defined as $\epsilon = \dfrac{T_1}{T_N * N}$ where $N$ is the number 
of MPI processes and $T_j$ is the wall time needed to complete the calculation with $j$ MPI processes. 
For a perfect scaling implementation $\epsilon$ is equal to one.
The parallel speedup with N processors is given by $S = T_N / T_1$.
Note that ``autoparal = 1`` will automatically change your ``job.sh`` script as well as the input file 
so that we can run the job in parallel with the optimal configuration required by the user. 
For example, you can use ``paral_kgb = 1`` in GS calculations and AbiPy will automatically set the values 
of ``npband``, ``npfft``, ``npkpt`` ... for you! 
Note that if no configuration fulfills the given condition, AbiPy will use the optimal configuration 
that leads to the highest parallel speedup (not necessarily the most efficient one).

``policy`` 
    This section governs the automatic parallelization of the run: in this case AbiPy will use 
    the ``autoparal`` capabilities of Abinit to determine an optimal configuration with 
    **maximum** ``max_ncpus`` MPI nodes. Setting ``autoparal`` to 0 disables the automatic parallelization. 
    Other values of autoparal are not supported*
