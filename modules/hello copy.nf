process hello {
    input:
    path input_file

    output:
    path 'hello_out.txt'

    script:
    """
    bash hello.sh ${input_file} > "hello_out.txt"
    """
}