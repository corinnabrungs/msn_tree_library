# metadata cleanup

```bash
pip install requirements
```

## run with prefect 2

```bash
prefect server start
```

Deployment:

- Start serving the flow locally by running metadata_cleanup_prefect.py, the service will then run
  and wait for jobs to be submitted
- or by running deployment in the terminal

```bash
prefect deploy metadata_cleanup_prefect:local-deploy
```

## Create and run a worker pool

- create worker pool with the name defined in the deployment (e.g., see metadata_cleanup_prefect.py
  main).

```bash
prefect work-pool create --type process local-work
prefect work-pool update --concurrency-limit 5 local-work
```

### start worker in pool to process
```bash
prefect work-pool update --concurrency-limit 5 local-work
prefect worker start --pool local-work
```

## Run jobs

Define jobs in jobs.py and run on prefect deployment.

