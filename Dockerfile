FROM jupyter/base-notebook

USER $NB_UID
RUN pip install ecell4
