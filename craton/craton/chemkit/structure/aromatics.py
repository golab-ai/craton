from copy import deepcopy

"""
    NA: total number of atoms of ring
    NDBI: the number of double bonds inside ring
    NDBE: the number of double bonds outside ring but connected with ring
    NO: the number of oxygen or sulfur atoms in ring
    NN3: the number of nitrogen atoms in ring with three connectivity 
    NN3R: the number of nitrogen atoms in multi ring with three connectivity
    NCP: the number of CH2 in ring connect with atom which has a double bond in ring
    NDBR: requested number of double bonds to forming planar ring
    Npie: the electron number contributed to forming large pi bond

    Na = atom_num
    Ndbi = double_bond_in_ring_num
    Ndbe = double_bond_out_ring_num
    No = o_or_s_number
    Nn3 = n_num_with_three_conn
    NN3R = n_num_with_three_conn_multi_ring
    NCP = n_num_CH2_connect_with_double_bond
    Ndbr = Na/2
    Npie=(Ndbi+No+Nn3)*2 #large_pi_bond

    aromatics:
        ar1: the standard aromatic ring, in which there are conjugate double-single bonds with out any Ndbe, No, Nn3
             and the number of conjugate electrons meet the 4n+2 rule, such as c1ccccc1(benzene), c1cnccc1(pyridine)
        ar2: the other standard aromatic ring, in which the ring member is 5, and there is one No or Nn3, without Ndbe
             such as o1cccc1(furan), [nH]1cccc1(pyrrole), [nH]cncc1(imidazole)
        ar3: plenary ring, in which there are No, Nn3, but without Ndbe
        ar4: plenary ring, in which there are Ndbe, with or without No, Nn3
        ar5: plenary ring, in which there are NCP
        ar6: the ring, in which there are conjugate double-single bonds, but the number of conjugate electrons not meet the 4n+2 rule
        #ar7: the ring, in which there are double bonds
        #ar8: the ring, in which there are Ndbe
        nonar:
"""

class Aromatics:
    """ """

    def __init__(self, molecule_obj, arr, total_ring_atoms):
        """
        m -> obj, Molecule object
        arr -> arr in ring_arr
              e.g.: ring_arr[0] -> [0,1,2,3,4,5,0]
        """
        self.m = deepcopy(molecule_obj)
        self.arr = arr
        self.total_ring_atoms = total_ring_atoms

    def get_total_atom_in_ring(self):
        """
        arr -> [5,4,3,2,1,0,5]

        return: 6
        """

        return len(self.arr[:-1])

    def get_double_bond_in_ring(self):
        """
        arr -> [0,1,2,3,4,5,0]

        e.g.:
        index: 0
        connect: [1,5,6]
        bond_type: ['2','1','1']

        return: 3
        """

        count = 0
        pair = []
        for index in self.arr[:-1]:
            atom = self.m.Atoms[index]
            connect = atom.connect
            #if hasattr(atom, "bond_type_old"):
            #    bond_type = atom.bond_type_old
            #else:
            #    bond_type = atom.bond_type
            if "2" not in atom.bond_type:
                continue
            else:
                for idx, item in enumerate(connect):
                    if item in self.arr[:-1] and atom.bond_type[idx] == "2":
                        if (index, item) not in pair and (item, index) not in pair:
                            pair.append((index, item))
                            count += 1
        return count

    def get_double_bond_out_ring(self):
        """
        arr -> [0,1,2,3,4,5,0]

        e.g.:
        index: 0
        connect: [1,5,6]
        bond_type: ['2','1','1']

        return: 0

        """

        count1 = 0
        count2 = 0
        pair = []
        for index in self.arr[:-1]:
            atom = self.m.Atoms[index]
            connect = atom.connect
            #if hasattr(atom, "bond_type_old"):
            #    bond_type = atom.bond_type_old
            #else:
            #    bond_type = atom.bond_type
            if "2" not in atom.bond_type:
                continue
            else:
                for idx, item in enumerate(connect):
                    if item not in self.arr[:-1] and atom.bond_type[idx] == "2":
                        if (index, item) not in pair and (item, index) not in pair:
                            if item in self.total_ring_atoms:
                                count2 += 1
                            else:
                                count1 += 1
                            pair.append((index, item))

        return count1, count2

    def double_bond_in_connect_atom(self,ii,bond_in_ring=False,atom_in_ring=True):
        if atom_in_ring:
            tmp = [jj for jj in self.m.Atoms[ii].connect if jj in self.arr[:-1]]
        else:
            tmp = self.m.Atoms[ii].connect
        tmp_tmp = []
        for jj in tmp:
            if bond_in_ring:
                tmp_tmp.extend([1 for nn,kk in enumerate(self.m.Atoms[jj].connect) 
                            if kk != ii and kk in self.arr[:-1] and self.m.Atoms[jj].bond_type[nn] in ["2","ar"]])
            else:
                tmp_tmp.extend([1 for nn,kk in enumerate(self.m.Atoms[jj].connect) 
                            if kk != ii and self.m.Atoms[jj].bond_type[nn] in ["2","ar"]])
            #if len(tmp_tmp) > 0:
            #    return True
        #return False
        return len(tmp_tmp)

    def get_special_atom_num_in_ring(self):
        """
        arr -> [0,1,2,3,4,5,0]

        """
        O_or_S = 0
        CH2 = 0
        NN3 = 0
        NN3_multi = 0
        for ii in self.arr[:-1]:
            elem = self.m.Atoms[ii].elem

            connect = len(self.m.Atoms[ii].connect)
            formal_charge = 0 if not hasattr(self.m.Atoms[ii], "formal_charge") else self.m.Atoms[ii].formal_charge


            if elem in ["O","S"]:
                flag = self.double_bond_in_connect_atom(ii)
                if flag >= 1:
                    O_or_S += 1
            elif elem == "C" and connect == 4:
                flag = self.double_bond_in_connect_atom(ii,bond_in_ring=True)
                if flag >= 2:
                    CH2 += 1

            if elem == "N" and connect - formal_charge == 3:
                ring_number = self.total_ring_atoms.count(ii)
                #ring_number = 1 # TO delete
                if ring_number > 1:
                    NN3_multi += 1
                else:
                    flag = self.double_bond_in_connect_atom(ii)
                    if flag >= 1:
                        NN3 += 1

        return O_or_S, CH2, NN3, NN3_multi


    def old_run_get_aromatics(
        self,
        atom_num,
        double_bond_in_ring_num,
        double_bond_out_ring_num,
        double_bond_out_ring_in_other_ring_num,
        o_or_s_number,
        n_num_with_three_conn,
        n_num_with_three_conn_multi_ring,
        n_num_CH2_connect_with_double_bond,
    ):
        planar_ring = atom_num / 2
        double_bond_num = double_bond_in_ring_num + 0.5 * (
            double_bond_out_ring_in_other_ring_num
            + double_bond_out_ring_num
            + o_or_s_number
            + n_num_with_three_conn
            + n_num_with_three_conn_multi_ring
        )
        #+ n_num_CH2_connect_with_double_bond
        large_pi_bond = (
            double_bond_in_ring_num + o_or_s_number + n_num_with_three_conn
        ) * 2 + double_bond_out_ring_in_other_ring_num

        if len(self.arr[:-1]) == 5:
            large_pi_bond += 2 * n_num_with_three_conn_multi_ring
            n_num_with_three_conn += n_num_with_three_conn_multi_ring
        else:
            large_pi_bond += n_num_with_three_conn_multi_ring

        if double_bond_num == planar_ring:
            if 0 == double_bond_out_ring_num:
                # AR1 or AR2 or AR4

                if 0 == (large_pi_bond - 2) % 4:
                    # AR1 or AR2(AR)
                    if 0 == o_or_s_number and 0 == n_num_with_three_conn:
                        return "ar1"
                    else:
                        return "ar2"
                else:
                    return "ar4"
            else:
                return "ar3"
                # AR3 or AR4
                # if 0 == o_or_s_number and 0 == n_num_with_three_conn:
                #    return "ar3"
                # else:
                #    return "ar4"
        else:
            if large_pi_bond >= planar_ring * 2 - 1 and large_pi_bond <= planar_ring * 2 + 1:
                return "ar5"
            else:
                if 0 == double_bond_out_ring_num:
                    if 0 == double_bond_in_ring_num:
                        return "nonar"
                    else:
                        return "ar8"
                else:
                    if 0 == double_bond_in_ring_num:
                        return "ar7"
                    else:
                        return "ar6"

    def run_get_aromatics(
        self,
        atom_num,
        double_bond_in_ring_num,
        double_bond_out_ring_num,
        double_bond_out_ring_in_other_ring_num,
        o_or_s_number,
        n_num_with_three_conn,
        n_num_with_three_conn_multi_ring,
        n_num_CH2_connect_with_double_bond,
    ):


        planar_ring = atom_num / 2
        double_bond_num = double_bond_in_ring_num + 0.5 * (
            double_bond_out_ring_in_other_ring_num
            + double_bond_out_ring_num
            + o_or_s_number
            + n_num_with_three_conn
            + n_num_with_three_conn_multi_ring
            + n_num_CH2_connect_with_double_bond
        )
        large_pi_bond = (
            double_bond_in_ring_num + o_or_s_number + n_num_with_three_conn
        ) * 2 + double_bond_out_ring_in_other_ring_num
        if len(self.arr[:-1]) == 5:
            large_pi_bond += 2 * n_num_with_three_conn_multi_ring
            n_num_with_three_conn += n_num_with_three_conn_multi_ring
        else:
            large_pi_bond += n_num_with_three_conn_multi_ring

        if double_bond_num == planar_ring:
            if planar_ring > 3:
                return "ar6"
            else:
                if 0 == n_num_CH2_connect_with_double_bond:
                    # ar1 ar2 ar3 or ar4
                    if 0 == double_bond_out_ring_num:
                        # ar1, ar2 or ar3
                        if planar_ring == 2.5:
                            if o_or_s_number + n_num_with_three_conn <= 1:
                                return "ar2"
                            else:
                                return "ar3"
                        else:
                            if 0 == o_or_s_number and 0 == n_num_with_three_conn:
                                return "ar1"
                            else:
                                return "ar3"
                    else:
                        return "ar4"
                else:
                    return "ar5"
        else:
            return "nonar"
            if 0 == double_bond_out_ring_num:
            #    if 0 == double_bond_in_ring_num:
                return "nonar"
            #    else:
            #        return "ar7"
            #else:
            #    return "ar8"


    def get_aromatics(self):
        atom_num = self.get_total_atom_in_ring()
        double_bond_in_ring_num = self.get_double_bond_in_ring()
        double_bond_out_ring_num, double_bond_out_ring_in_other_ring_num = self.get_double_bond_out_ring()
        o_or_s_number, n_num_CH2_connect_with_double_bond,n_num_with_three_conn,n_num_with_three_conn_multi_ring = self.get_special_atom_num_in_ring()
        #n_num_with_three_conn, n_num_with_three_conn_multi_ring = self.get_three_bond_with_N_in_ring()
        return self.run_get_aromatics(
            atom_num,
            double_bond_in_ring_num,
            double_bond_out_ring_num,
            double_bond_out_ring_in_other_ring_num,
            o_or_s_number,
            n_num_with_three_conn,
            n_num_with_three_conn_multi_ring,
            n_num_CH2_connect_with_double_bond,
        )

    def __call__(self):
        func = self.get_aromatics
        func()
