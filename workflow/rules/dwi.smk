def get_import_dwi_cmd(wildcards,input,output):
    bval = sorted(glob( input.nii_folder + f'/*cfmmDtiEpi_10B0_30B1k_60B2k_Marm_ISO400.bval'))[-1] #get last one only
    prefix=bval[:-5]
    cmds=[]
    nii_no_gz = output.nii[:-3]

    cmds.append(f'cp {prefix}.nii {nii_no_gz}')
    cmds.append(f'gzip {nii_no_gz}')
    cmds.append(f'cp {prefix}.bval {output.bval}')
    cmds.append(f'cp {prefix}.mvec {output.bvec}')
    cmds.append(f'cp {prefix}.json {output.json}')
    return ' && '.join(cmds)



rule import_dwi_files:
    input:
        nii_folder='niftis/sub-{subject}',
    params:
        cmd=get_import_dwi_cmd
    output:
        nii=bids(root='bids',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.nii.gz'),
        bval=bids(root='bids',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.bval'),
        bvec=bids(root='bids',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.bvec'),
        json=bids(root='bids',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.json')
    shell:
        "{params.cmd}"
     

def get_import_dwi_rev_cmd(wildcards,input,output):
    nii = sorted(glob( input.nii_folder + f'/*cfmmDtiEpi_5B0_RVPhase_Marm_ISO400.nii'))[-1] #get last one only
    prefix=nii[:-4]
    cmds=[]
    nii_no_gz = output.nii[:-3]

    cmds.append(f'cp {prefix}.nii {nii_no_gz}')
    cmds.append(f'gzip {nii_no_gz}')
    cmds.append(f'cp {prefix}.json {output.json}')

    #need to create bval and bvec file
    cmds.append(f'echo "0" > {output.bval}')
    cmds.append(f'echo "0\n0\n0\n" > {output.bvec}')


    return ' && '.join(cmds)

   
rule import_dwi_rev_files:
    input:
        nii_folder='niftis/sub-{subject}',
    params:
        cmd=get_import_dwi_rev_cmd
    output:
        nii=bids(root='bids',
                subject='{subject}',
                acq='revb0',
                datatype='dwi',
                suffix='dwi.nii.gz'),
        bval=bids(root='bids',
                subject='{subject}',
                acq='revb0',
                datatype='dwi',
                suffix='dwi.bval'),
        bvec=bids(root='bids',
                subject='{subject}',
                acq='revb0',
                datatype='dwi',
                suffix='dwi.bvec'),
        json=bids(root='bids',
                subject='{subject}',
                acq='revb0',
                datatype='dwi',
                suffix='dwi.json')
    shell:
        "{params.cmd}"
        

