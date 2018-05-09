Hosted at Github
================

https://github.com/tgen/jetstream

Contributing
-------------

All contributions are welcome.
Comments, suggeestions, pull-requests, and issues can be submitted on `Github`_.

Development install
-------------------

Best practices for development are to clone the source repo, create a
virtual environment and install as an editable link with Pip. This will
cause (most) changes made in the source code to be testable without
resintalling. Here is an example command list:

::

   git clone https://github.com/tgen/jetstream.git
   virtualenv jetstream_venv
   source jetstream_venv/bin/activate
   pip3 install -e jetstream/

.. _Github: https://github.com/tgen/jetstream/issues