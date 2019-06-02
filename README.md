# ecell4

This provides pure Python libraries for E-Cell System version 4.

## Quick Start (with Docker)

```
docker run -p 8888:8888 ecell/ecell4
```
This command pulls the ecell/ecell4 image from Docker Hub.

It then starts a container running a Jupyter Notebook server and exposes the server on host port 8888.
(ecell4 is installed in this server environment.)

The server logs appear in the terminal.

Visiting `http://<hostname>:8888/?token=<token>` in a browser loads the Jupyter Notebook dashboard page, where hostname is the name of the computer running docker and token is the secret token printed in the console.

## Installation

### Any Linux (almost)

```
pip3 install ecell4
```

### Windows, Mac

Please refer to [ecell4-base INSTALL.md](https://github.com/ecell/ecell4-base/blob/master/INSTALL.md)
