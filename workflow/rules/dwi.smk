def get_import_dwi_cmd(wildcards,input,output):
    rev_nii = sorted(glob( input.nii_folder + f'/*cfmmDtiEpi_5B0_RVPhase_Marm_ISO400.nii'))[-1] #get last series only
    rev_prefix=rev_nii[:-4]
    cmds=[]

    cmds.append(f'fslsplit {rev_prefix}.nii split_rev_dwi -t')
    cmds.append(f'im=`ls split_rev_dwi*nii.gz | tail -n 1`')
    cmds.append(f'c3d $im -orient RAI -o {output.rev_nii}')  #take last b0 image only
    cmds.append(f'cp {rev_prefix}.json {output.rev_json}')

    #header affines are different for dti and rev_dti (arrrghhhh!) - use fslswapdim first

    nii = sorted(glob( input.nii_folder + f'/*cfmmDtiEpi_10B0_30B1k_60B2k_Marm_ISO400.nii'))[-1] #get last one only
    prefix=nii[:-4]
    cmds.append(f'fslsplit {prefix}.nii split_dwi -t')
    cmds.append(f'for im in `ls split_dwi*nii.gz`; do fslswapdim $im -x -y z $im; c3d $im -orient RAI -o $im; done')  
    cmds.append(f'fslmerge -t {output.nii} split_dwi*.nii.gz')
    cmds.append(f'cp {prefix}.json {output.json}')


    #need to create bval and bvec file for rev ph enc
    cmds.append(f'echo "0" > {output.rev_bval}')
    cmds.append(f'echo "0\n0\n0\n" > {output.rev_bvec}')





    return ' && '.join(cmds)



rule import_dwi_files:
    """ except for bval and bvec (that comes from bruker private tags)"""
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
        json=bids(root='work',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.json'),
        rev_nii=bids(root='bids',
                subject='{subject}',
                acq='revb0',
                datatype='dwi',
                suffix='dwi.nii.gz'),
        rev_bval=bids(root='bids',
                subject='{subject}',
                acq='revb0',
                datatype='dwi',
                suffix='dwi.bval'),
        rev_bvec=bids(root='bids',
                subject='{subject}',
                acq='revb0',
                datatype='dwi',
                suffix='dwi.bvec'),
        rev_json=bids(root='work',
                subject='{subject}',
                acq='revb0',
                datatype='dwi',
                suffix='dwi.json')
    shadow: 'minimal'
    shell:
        "{params.cmd}"
     
rule extract_bruker_info_dwi:
    input:
        nii_folder='niftis/sub-{subject}',
        dcm_folder='dicoms/sub-{subject}'
    output:
        bval=bids(root='work',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.bval'),
        bvec=bids(root='work',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.bvec'),
        bmat=bids(root='work',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.bmat'),
        method_json=bids(root='work',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.method.json'),
    script: '../scripts/extract_bruker_info_dwi.py'

rule dwi_gradient_fix:
    input:
        nii=bids(root='bids',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.nii.gz'),
        bval=bids(root='work',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.bval'),
        bvec=bids(root='work',
                subject='{subject}',
                acq='multishell',
                datatype='dwi',
                suffix='dwi.bvec'),
    output:
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
    shadow: 'minimal'
    shell: 'dwigradcheck {input.nii} -fslgrad {input.bvec} {input.bval} -export_grad_fsl {output.bvec} {output.bval}'


def get_phase_encoding_direction(wildcards):
    if wildcards.acq == 'multishell':
        return 'i'
    
    if wildcards.acq == 'revb0':
        return 'i-'

rule finalize_json_dwi:
    input:
        json=bids(root='work',
                subject='{subject}',
                acq='{acq}',
                datatype='dwi',
                suffix='dwi.json'),
    params:
        phase_enc = get_phase_encoding_direction
    output:
        json=bids(root='bids',
                subject='{subject}',
                acq='{acq}',
                datatype='dwi',
                suffix='dwi.json'),
    run:
        import json
        with open(input.json, 'r') as f:
            data = json.load(f) 

        data['PhaseEncodingDirection']=params.phase_enc

        with open(output.json, 'w') as f:
            json.dump(data, f, indent=4)



