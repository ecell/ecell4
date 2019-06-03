# ecell4

This provides pure Python libraries (DB integration and miscellaneous utility functions) for [E-Cell System version 4](https://github.com/ecell/ecell4-base).

Installation
------------

### Windows

- Install Miniconda with **Python3.7 64-bit** from http://conda.pydata.org/miniconda.html
- Run the following commands in the command prompt.

    ```shell
    conda update conda python pip
    conda uninstall hdf5
    conda clean -a
    conda install hdf5 notebook
    pip install ecell4
    ```

If you need to create movie with E-Cell4, please install [ffmpeg windows build](http://ffmpeg.zeranoe.com/builds/) and add its path to your **USER** PATH enviromental variable.

### Mac

- Install Miniconda with **Python3.7 64-bit** from http://conda.pydata.org/miniconda.html
- Run the following commands in the Terminal app.

    ```shell
    ~/miniconda3/bin/conda update conda python pip
    ~/miniconda3/bin/conda install notebook
    ~/miniconda3/bin/pip install ecell4
    # If you need to create movie, install ffmpeg with homebrew
    brew install ffmpeg
    ```

### Any Linux (almost)

You do not need to install Miniconda on Linux.

```shell
python3 -m pip install ecell4
# If you need to create movie[, and use Ubuntu]
apt install ffmpeg
```



### Docker

```
docker run -p 8888:8888 ecell/ecell4
```
This command pulls the ecell/ecell4 image from Docker Hub.

It then starts a container running a Jupyter Notebook server and exposes the server on host port 8888.
(ecell4 is installed in this server environment.)

The server logs appear in the terminal.

Visiting `http://<hostname>:8888/?token=<token>` in a browser loads the Jupyter Notebook dashboard page, where hostname is the name of the computer running docker and token is the secret token printed in the console.
