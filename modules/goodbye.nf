process goodbye {
    input:
    path input_file

    output:
    path 'goodbye_out.txt'

    script:
    """
    bash goodbye.sh ${input_file} > "goodbye_out.txt"
    """
}