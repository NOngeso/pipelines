#!/usr/bin/python

import web
import kang

file_fa = 'Creinhardtii_281_v5.0.fa'

dicHD2seq = kang.Fasta2dic(file_fa)



urls = ('/', 'index')
render = web.template.render('templates/')

def get_chopped(strin):
    window = 100
    out_list = []
    for i in range(len(strin)/window+1):
        try:
            out_list.append(strin[(i)*50:(i+1)*50])
        except IndexError:
            out_list.append(strin[(i)*50:])
    return out_list

class index:
    form = web.form.Form(
                 web.form.Textbox('chromosome',web.form.notnull, description="Chromosome number"),
                 web.form.Textbox('left',web.form.notnull, description="Left postion"),
                 web.form.Textbox('right',web.form.notnull, description="Right postion"),
                 web.form.Button('Submit'),
                )
    def GET(self):
        #i = web.input(name=None)
        chromosome_list = set(dicHD2seq.keys())
        
        form = self.form()
        return render.index(chromosome_list,form)
    def POST(self):
        form = self.form()
        if not form.validates():
            return render.index(form)
        chromosome = form.d.chromosome.strip()
        left = int(form.d.left)-1
        right = int(form.d.right)
        seq = dicHD2seq[chromosome][left:right]
        seqs = get_chopped(seq)
        return render.showseq(seqs, form)

if __name__ == "__main__": 
    app = web.application(urls, globals())
    app.run()    

