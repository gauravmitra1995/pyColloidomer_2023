def define_interactions(table,nl,softrepulsion,type_list,softV_epsilon_dict,softV_rcut_dict,r_dict):

    for x in range(len(type_list)):
        for y in range(len(type_list)):
            if(x<=y):
                first=str(type_list[x])
                second=str(type_list[y])
                key=first+','+second
                if key in softV_epsilon_dict.keys():
                   table.pair_coeff.set(first,second,func=softrepulsion,rmin=0,rmax=1.05*softV_rcut_dict[key],\
                   coeff=dict(epsilon=softV_epsilon_dict[key],ron=0.1*softV_rcut_dict[key],rcut=softV_rcut_dict[key]))
                else:
                    rcut=2*r_dict[str(type_list[2])]
                    table.pair_coeff.set(first,second,func=softrepulsion,rmin=0,rmax=1.05*rcut,\
                    coeff=dict(epsilon=0,ron=-0.1*rcut,rcut=-rcut))
