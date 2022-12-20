import os.path
import os


list = os.listdir('/home/ytirlet/Documents/yael/prokka_pub/')
dico_prefix = {'Bifidobacterium':'B','Lactobacillus':'Lb', 'L':'Lco', 'Lactococcus':'Lco',  'Leuconostoc':'Leu', 'Pediococcus':'Pe', 'Propionibacterium':'Pr', 'P':'Pr', 'Streptococcus':'St'}
dico_suffix = {'Complete':'C', 'complet':'C', 'C':'C', 'Scaffold':'S', 'S':'S', 'Plasmid':'P', 'plasmide':'P', 'P':'P'}
for name in list :
    parts = name.split('_')
    prefix = parts[0]
    middle = parts[1:-1]
    suffix = parts[-1]
    if prefix in dico_prefix :
        new_name = dico_prefix[prefix] + '_'
    else :
        new_name = prefix + '_'
    for i in range(len(middle)) :
        new_name += middle[i] + '_'
    if suffix in dico_suffix :
        new_name += dico_suffix[suffix]
    else :
        new_name += suffix + '_S'
    if '.' in new_name or ':' in new_name :
        while '.' in new_name or ':' in new_name :
            new_name = new_name.replace('.','-')
            new_name = new_name.replace(':','-')
    os.system("mv -fv /home/ytirlet/Documents/yael/prokka_pub/"+name+"/*.gbk /home/ytirlet/Documents/yael/prokka_pub/"+name+"/"+new_name+".gbk")
    os.system("mv -fv /home/ytirlet/Documents/yael/prokka_pub/"+name+"/ /home/ytirlet/Documents/yael/prokka_pub/"+new_name+"/")

