import os
import argparse
import subprocess
import threading
import uuid
import socket

import concurrent.futures

worker_script = """#!/bin/bash

slog=/tmp/worker-startup.log
workerdownload=https://go.parallel.works/fixed/lib/worker.pl

echo "RUNNING STARTUP SCRIPT"

slog()
{{
  TS=$(date +%Y.%m%d.%H%M%S).$(date +%N | sed -e 's/......$//')
  echo $TS $* >> $slog
}}

slog gcloud Linux worker startup script log

# --- Handling threads --- #
disable_ht () {{
    for cpunum in $(cat /sys/devices/system/cpu/cpu*/topology/thread_siblings_list | cut -s -d, -f2- | tr ',' '\n' | sort -un)
    do
	echo 0 > /sys/devices/system/cpu/cpu$cpunum/online
    done
}}

#if __DISABLE_HT__; then
#    disable_ht
#fi

id -u pwuser
userexists=$(echo $?)

mount_nfs () {{
    sudo mkdir -p ${{nfs}}
    #storageacct="$(echo ${{pwuser}} | cut -c1-4)$(echo ${{pool}} | cut -c1-4)$(echo ${{uuid}} | sed 's/-//g')strga"
    #storagekey=$(az storage account keys list --resource-group "${{group}}" --account-name ${{storageacct}} --query "[0].value" | tr -d '"')
    storageacct=__STORAGEACCT__
    storagekey=__STORAGE_KEY__
    stshname=__STSHNAME__
    sudo mount -t cifs //${{storageacct}}.file.core.windows.net/$(echo ${{stshname}}) $nfs -o vers=3.0,username=${{storageacct}},password=${{storagekey}},dir_mode=0777,file_mode=0777,serverino
    exit_status=$?
    until [ ${{exit_status}} -eq 0 ]; do
	sleep ${{group_info_freq}}
	#storagekey=$(az storage account keys list --resource-group "${{group}}" --account-name ${{storageacct}} --query "[0].value" | tr -d '"')
	sudo mount -t cifs //${{storageacct}}.file.core.windows.net/$(echo ${{groupname}})-fileshare $nfs -o vers=3.0,username=${{storageacct}},password=${{storagekey}},dir_mode=0777,file_mode=0777,serverino
	exit_status=$?
    done
    sudo chmod 0777 $nfs
}}

# Using az client
az_create_hostfile () {{
    groupid=/subscriptions/d53028c1-2dd1-47c8-bb20-8eb1f5ea9338/resourceGroups/${{group}}/providers/Microsoft.Compute/virtualMachines/
    vmsids=$(seq -f "${{groupid}}${{groupname}}-%02g-$(printf %02d ${{groupsize}})-${{groupzone}}" 1 ${{groupsize}})
    az vm list-ip-addresses --ids ${{vmsids}} --query [].virtualMachine.network.publicIpAddresses[].ipAddress -o tsv > ${{hostfile}}
    
}}

# Using naming convention and curl
create_hostfile () {{
    rm -f ${{hostfile}}
    vms=$(seq -f "${{groupname}}-%02g-$(printf %02d ${{groupsize}})-${{groupzone}}" 1 ${{groupsize}})
    for vm in ${{vms}}; do
	echo "Adding IP for $vm"
	until [ -f ${{nfs}}/${{vm}}.ip ]; do
	    sleep ${{group_info_freq}}
	done
	cat ${{nfs}}/${{vm}}.ip
	cat ${{nfs}}/${{vm}}.ip >> ${{hostfile}}
        vmip=$(cat ${{nfs}}/${{vm}}.ip)	
	ssh-keyscan ${{vmip}} >> ~/.ssh/known_hosts
	ssh -o "StrictHostKeyChecking no" -t ${{vmip}} "echo Hello from \$(hostname)"
    done
}}
   
generate_sshkeys () {{
    if [[ "${{nodenumber}}" == "01" ]]; then
        ssh-keygen -f ~/.ssh/id_rsa -t rsa -N ''
        cp ~/.ssh/id_rsa.pub ${{nfs}}
    fi
    until [ -f ${{nfs}}/id_rsa.pub ]; do
	sleep ${{group_info_freq}}
    done
    pubkey=$(cat ${{nfs}}/id_rsa.pub)
    echo $pubkey >> ~/.ssh/authorized_keys
}}

	
if [[ "$userexists" == "1" ]];then

echo "Creating pwuser login..."
slog "Creating pwuser login..."

# Create Linux User and cd to its home directory

sudo mkdir -p /pworks
sudo addgroup --gid 9001 pwuser

sudo useradd -d /pworks -M -N -s /bin/bash -c "Parallel Works app-run user" -g pwuser -G docker -u 9001 pwuser

cd /pworks
# mkdir -p workerlog # Not needed?  Put workerlogs in $PWD of /pworks/u  => or just /pworks for now

# Set Docker PW credentials

#if [[ -d /pworks/.docker ]]; then
#    mv /pworks/.docker /pworks/.docker.save
#fi
#if [[ -s /root/.docker/config.json ]]; then
#  cp -rp /root/.docker /pworks/.docker
#fi    

sudo chown -R pwuser /pworks
sudo chgrp -R pwuser /pworks

WORKER_LOGGING_LEVEL=INFO;
printenv >>$slog

#cat <<END >>/etc/sudoers 
#pwuser ALL=(root) NOPASSWD:ALL
#END

:<<ENDALL
# THERE IS A PROBLEM WITH THIS - CANNOT BUILD IMAGES FROM THIS - FSTAB DEFINITIONS MESS UP FUTURE RESTARTS

# Enable pwuser to mount and format disks
cat <<END >>/etc/sudoers                 # FIXME: May not be needed
pwuser ALL=(root) NOPASSWD:/bin/mount
pwuser ALL=(root) NOPASSWD:/bin/umount
pwuser ALL=(root) NOPASSWD:/sbin/mke2fs
pwuser ALL=(root) NOPASSWD:/sbin/mkfs.ext4
END

# Enable pwuser to format (mkfs) and mount file systems on /dev/sdb and beyond

# chown pwuser /dev/sd{{b,c,d,e,f,g}} # FIXME: these dont exist at boot. Use sudo mkfs.ext4
# chgrp pwuser /dev/sd{{b,c,d,e,f,g}}

mkdir -p /pworks/mnt/{{b,c,d,e,f,g}}
chown pwuser /pworks/mnt/{{b,c,d,e,f,g}}
chgrp pwuser /pworks/mnt/{{b,c,d,e,f,g}}

cat <<END >>/etc/fstab
/dev/sdb /pworks/mnt/b ext4 discard,user,exec,defaults 0 2
/dev/sdc /pworks/mnt/c ext4 discard,user,exec,defaults 0 2
/dev/sdd /pworks/mnt/d ext4 discard,user,exec,defaults 0 2
/dev/sde /pworks/mnt/e ext4 discard,user,exec,defaults 0 2
/dev/sdf /pworks/mnt/f ext4 discard,user,exec,defaults 0 2
/dev/sdg /pworks/mnt/g ext4 discard,user,exec,defaults 0 2
END
ENDALL

fi

start_worker () {{
    # Download worker.pl script and set dir/file permissions
    slog "fetching worker script: wget -O worker.pl $workerdownload"
    wget -O worker.pl $workerdownload  >>$slog
    chmod +x worker.pl
    ls -l worker.pl >>$slog
    
    # Launch worker as userid pwuser

    worker_cmd="./worker.pl -l INFO {service_url} {worker_name} $PWD"
    slog "launching worker command: $worker_cmd"

    #sudo -E -u pwuser bash -c "$worker_cmd >& worker.out ; echo worker.pl exit code: $? >>worker.out" &
    # causing issues with many workflow so temporarily running as user worker starts on
    $worker_cmd >& worker.out ; echo worker.pl exit code: $? >>worker.out &

    # FIXME: Need to dissociate worker from startup script? setsid etc?
    # FIXME: Copy worker logs to GCS on worker exit?  Always?  Conditionally?

    # For reference:
    #
    #       worker.pl [-l <loglevel>] [-w <maxwalltime>] [-c <concurrency>] <serviceURL> <blockID> <logdir>
}}

multihost=false
if [[ "${{multihost}}" == "true" ]]; then
    echo "Bootstrapping Group"
    hostname=$(cat /proc/sys/kernel/hostname)
    nfs=/mnt/shared
    pwuser=$(echo ${{hostname}} | cut -d "-" -f2)
    pool=$(echo ${{hostname}} | cut -d "-" -f3)
    uuid=$(echo ${{hostname}} | cut -d "-" -f4-5)
    groupname=w-${{pwuser}}-${{pool}}-${{uuid}}
    nodenumber=$(echo ${{hostname}} | cut -d'-' -f6)
    groupsize=$(echo ${{hostname}} | cut -d'-' -f7)
    groupzone=$(echo ${{hostname}} | cut -d "-" -f8-)
    hostfile="${{nfs}}/${{hostname}}.hostfile"
    group_info_freq=20
    mount_nfs
    echo $(curl ifconfig.me) > ${{nfs}}/${{hostname}}.ip
    if [[ "${{nodenumber}}" == "01" ]]; then
        echo I am the lead node
	create_hostfile
    fi
    generate_sshkeys
    if [[ "${{nodenumber}}" == "01" ]]; then
	start_worker
    fi
else
    start_worker
fi"""


def start_cluster(serviceport=4050, localport=4051, logs_dir=None, swift_bin_dir=None, stats=False):
    logs_dir = logs_dir or '/tmp/cluster_logs/coasters'
    swift_bin_dir = (swift_bin_dir + os.sep) if swift_bin_dir else ''
    stats = '-stats' if stats else ''
    print(swift_bin_dir)
    subprocess.call(f'kill -9 $(netstat -tulpn | grep {serviceport} | awk \'{{print $7}}\' | '
                     'cut -d "/" -f1) > /dev/null 2>&1', shell=True)
    subprocess.call(f'kill -9 $(netstat -tulpn | grep {localport} | awk \'{{print $7}}\' | '
                     'cut -d "/" -f1) > /dev/null 2>&1', shell=True)
    subprocess.call(f'rm {os.path.dirname(logs_dir)} -R > /dev/null 2>&1', shell=True)
    subprocess.call(f'{swift_bin_dir}coaster-service -p {serviceport} -localport {localport} {stats} -nosec '
                    f'-passive -logdir {logs_dir}', shell=True)



def start_remote_worker(worker_fqd, service_url, worker_name=None):
    worker_name = worker_name or str(uuid.uuid4())
    worker_start_script = worker_script.format(service_url=service_url, worker_name=worker_name)
    worker_start_script_path = f'.{worker_name}.worker_start.bash'
    with open(worker_start_script_path, 'w') as script_out:
        script_out.write(worker_start_script)
    
    with open(os.devnull, 'w') as dev_null:
        subprocess.call(['rsync', '-avzP', worker_start_script_path, f'{worker_fqd}:~'], stdout=dev_null)
    subprocess.call(f'ssh {worker_fqd} "sudo bash {worker_start_script_path}" >/dev/null 2>&1 &', shell=True)
    
    try:
        os.remove(worker_start_script_path)
    except:
        pass  # Fail silently



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fqdns', nargs='*')
    parser.add_argument('--swift-bin-dir')
    parser.add_argument('--show-stats', action='store_true')


    # parser.add_argument('-u', '--username', required=True)
    # parser.add_argument('--resource-group')
    # parser.add_argument('--image')
    # parser.add_argument('--size')
    # parser.add_argument('--azure-cli-container', help='Path to singularity container with azure cli in path')
    args = vars(parser.parse_args())
    
    # with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
    #     executor.submit(start_cluster)
    
    t = threading.Thread(target=start_cluster, kwargs={'swift_bin_dir': args['swift_bin_dir'], 'stats': args['show_stats']})
    t.start()
    
    
    local_ip = socket.gethostbyname(socket.gethostname())
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(args['fqdns'])) as executor:
        for fqdn in args['fqdns']:
            executor.submit(start_remote_worker, fqdn, f'http://{local_ip}:4051')
            print(f'Started worker node at {fqdn}')
            # start_remote_worker(fqdn, f'http://{local_ip}:4051')
    
    # start_cluster(swift_bin_dir=args['swift_bin_dir'])
    # print('done')




