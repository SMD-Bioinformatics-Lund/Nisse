process scout_yaml {
    input:
    path input_file

    output:
    path 'hello_out.txt'

    script:
    """
    bash hello.sh ${input_file} > "hello_out.txt"
    """
}