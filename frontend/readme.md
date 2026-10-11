# TBRC Frontend

This directory contains the Shiny frontend for TBRC.

It is responsible for:

- collecting run metadata
- selecting barcode presets
- choosing processing options
- launching upload-mode or server-mode submissions
- showing submission progress

## Main Files

- [`ui.R`](frontend/ui.R)
- [`server.R`](frontend/server.R)
- [`global.R`](frontend/global.R)
- [`server_config.json`](frontend/server_config.json)
- [`scripts/`](frontend/scripts)

## Machine-Specific Configuration

Most machine-specific settings now live in [`server_config.json`](frontend/server_config.json), including:

- storage server raw-data paths
- storage server results paths
- storage server SSH ports
- transfer mode
- cluster hostnames
- pipeline root
- backend environment name
- IGBlast paths

To move the app to a new server, update the JSON config rather than editing the R code or shell scripts.

## Required R Packages

At minimum:

```r
install.packages(c("shiny", "DT", "jsonlite", "shinyalert", "shinyWidgets"))
```

The app has a fallback path when `shinyWidgets` is missing, but installing it is recommended.

## Local Run

```bash
cd /path/to/TBRC/frontend
Rscript -e "shiny::runApp('.', host='127.0.0.1', port=8080)"
```

## Launcher Scripts

The launcher scripts in [`scripts/`](frontend/scripts) resolve paths relative to the repository, which makes the frontend easier to move between servers.

Important behavior:

- server-mode submissions validate the requested storage folder before launching the backend
- upload-mode submissions write files under the configured backend input area
- archive delivery is controlled by backend settings carried through `server_config.json`
