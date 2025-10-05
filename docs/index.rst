Home
====

.. Hell world

.. raw:: html

    <div style="text-align:center; padding-top: 1rem">
        <img src="_static/backpop.png" alt="BackPop Logo" style='width:100%; max-width: 700px'>
        <h4>A tool to sample the joint distributions of initial binary parameters and binary interaction assumptions</h4>
    </div>


.. test some code highlighting

.. code-block:: python

    from backpop import BackPop
    bp = BackPop()
    bp.sample(1000)
    bp.plot_distributions()

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :hidden:

   pages/install
   pages/getting_started
   pages/tutorials
   pages/modules
   pages/cite
