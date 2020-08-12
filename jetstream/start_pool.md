## Starting a pool of workers with ``start_pool.py``

### Pre-requisites
Before starting a pool, the user must first spin up the desired number of worker instances as well as have ``ssh`` access to all worker instances. It would be ideal for the connections to already be in the ``known_hosts`` files to make the process completely non-interactive.

### Starting the pool
From the desired client node, run ``start_pool.py`` to start the pool (in a GNU ``screen`` or as a background process if this is for a long running pipeline). The program has the following options:

* ``fqdns`` (positional) A space separated list of ``<username>@<hostname>`` strings, representing the necessary information to ``ssh`` into each worker node. There should be one list entry for each worker that is to be included in the pool.
* ``--swift-bin-dir`` The path to the directory containing the ``coaster-service`` swift binary. If not provided, ``coaster-service`` is assumed to be in path.
* ``--show-stats`` If provided, the program will show relevant pool statistics.

For example, given a pool of four workers:
``python start_pool.py --show-stats --swift-bin-dir /path/to/swift/bin user@1.0.0.1 user@1.0.0.2 user@1.0.0.3 user@1.0.0.4``

### Using the pool to run a Jetstream pipeline
When the pool is started, the program will present a _coaster service URL_:
* If the ``--show-stats`` option is used, the top-most section will display a line labeled _URL_
* If the ``--show-stats`` option is not used, the output to the terminal will contain a line ``Started coaster service:``

The _coaster service URL_ should be given to the Jetstream configuration file at ``~/.config/jetstream/config.yaml``. The relevant section will look like the following:

```
backends:
  cloud_swift:
    (): jetstream.backends.cloud.CloudSwiftBackend
    pw_pool_name: null
    pw_api_key: null
    azure_params:
      az_sif_path: '/home/user/azure-cli_latest.sif'
      az_storage_account_name: '<az_account_name>'
      az_storage_account_key: '<az_account_key>'
    pool_info:
      cpus_per_worker: 8
      workers: 2
      serviceurl: <coaster_service_URL>
```

If the user leaves ``backends.cloud_swift.pw_pool_name`` and ``backends.cloud_swift.pw_api_key`` are both set to ``null``, then the cloud backend will assume a pool has been manually created and attempt to use the definitions in ``backends.cloud_swift.pool_info``. 
* ``cpus_per_worker`` are the number of CPUs available in each worker (which is, for the moment, assumed to be homogeneous )
* ``workers`` is the number of workers in the pool
* ``serviceurl`` is the _coaster service URL_

Once all this configuration is in place, run the Jetstream pipeline with the ``cloud_swift`` backend and execution will begin on the workers in the pool.