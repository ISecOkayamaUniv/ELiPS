import os
ROOT_PATH = '/Users/hattoridaichi/workspace/elips_bitbucket/include/ELiPS'


def print_func_list(file_path):
    with open(file_path) as f:
      line_data=[]
      for s_line in f:
        if ('extern' in s_line) & ('(' in s_line):
          line_data.append(file_path.replace(ROOT_PATH + '/','')+','+s_line.replace('extern ', '').replace('void ','').replace('int ','').replace('\n',''))
          # print(file_path.replace(ROOT_PATH + '/','')+','+s_line.replace('extern ', '').replace('void ','').replace('int ','').replace('\n',''))
      outdata=[]
      for i in range(0,len(line_data)) :
        outdata.append(line_data[i].split("(")[0])

      for i in range(0,len(outdata)):
        print(outdata[i])

def recursive_file_check(path):
    if os.path.isdir(path):
        # directoryだったら中のファイルに対して再帰的にこの関数を実行
        files = os.listdir(path)
        for file in files:
            recursive_file_check(path + "/" + file)
    else:
        # fileだったら処理
        print_func_list(path)


recursive_file_check(ROOT_PATH)




# path = '/Users/hattoridaichi/workspace/elips_bitbucket/include/ELiPS/fp.h'
# with open(path) as f:
#   for s_line in f:
#     if 'extern' in s_line:
#       print(s_line.replace('extern ', ''))
#       print()
