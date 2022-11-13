def write_file(header, content, file_name):
    #Used to write files based on the input content
    with open(file_name + '.csv' if file_name[-4:] != '.csv' else file_name, 'w') as csvfile:
        for el in header:
            csvfile.write(el) if el != header[0] else csvfile.write(f'{el}') if el[0] != ':' else csvfile.write(el)
            csvfile.write(',') if el != header[-1] else csvfile.write('\n')
            
        for el_list in content:
            for el in el_list: 
                csvfile.write(str(el))
                csvfile.write(',') if el != el_list[-1] else csvfile.write('\n')