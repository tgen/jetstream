from jetstream import plugins

def slurm(plugin):
    #cmd = "srun {}".format(plugins[plugin])
    cmd = "{}".format(plugins[plugin])
    return cmd
