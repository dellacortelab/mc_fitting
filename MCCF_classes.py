import numpy as np


def lagrangian_strain(dfgrd1):
    '''calculate lagrangian strain from deformation gradient'''
    strain = np.zeros((3, 3), float)
    delta = np.zeros((3, 3), float)
    delta[0, 0] = 1
    delta[1, 1] = 1
    delta[2, 2] = 1

    for i in range(3):
        for j in range(3):
            for k in range(3):
                strain[i, j] += dfgrd1[k, i] * dfgrd1[k, j]
            strain[i, j] = 1/2. * (strain[i, j] - delta[i, j])
    return strain


def cauchy_green_strain(dfgrd1):
    '''calculate right cauchy green strain from deformation gradient'''
    strain = np.zeros((3, 3), float)
    delta = np.zeros((3, 3), float)
    delta[0, 0] = 1
    delta[1, 1] = 1
    delta[2, 2] = 1

    for i in range(3):
        for j in range(3):
            for k in range(3):
                strain[i, j] += dfgrd1[k, i] * dfgrd1[k, j]
    return strain


def pk2(cauchy, dfgrd1):
    '''calculate the 2 Piola Kirchhoff Stress fom Cauchy Stress'''

    pk2 = np.zeros((3, 3), float)
    idfgrd1 = umat.invert(dfgrd1, 3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for r in range(3):
                    pk2[i, j] += idfgrd1[i, r] * cauchy[r, k] * idfgrd1[j, k]
    return pk2


def pk2_voigt(cauchy, dfgrd1):
    '''calculate the 2 Piola Kirchhoff Stress fom Cauchy Stress, return voigt'''

    pk2 = np.zeros((3, 3), float)
    idfgrd1 = umat.invert(dfgrd1, 3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for r in range(3):
                    pk2[i, j] += idfgrd1[i, r] * cauchy[r, k] * idfgrd1[j, k]
    return [pk2[0, 0], pk2[1, 1], pk2[2, 2], pk2[0, 1], pk2[0, 2], pk2[1, 2]]


def mat6n33a(a6):
    '''! c....transfer 6-dim. Voigt notation to 3x3 matrix'''
    b9 = np.zeros((3, 3), float)
    for i in range(3):
        b9[i, i] = a6[i]

    b9[0, 1] = a6[3]
    b9[1, 0] = a6[3]
    b9[1, 2] = a6[5]
    b9[2, 1] = a6[5]
    b9[2, 0] = a6[4]
    b9[0, 2] = a6[4]
    return b9


def load_experimental_data(filename):
    '''reads an x y text file'''
    x, y = np.loadtxt(filename, float)
    return x, y


class Material_model():
    '''Material model calss. It requieres and instance of the properties class
       to be initalized. It provides an interface to set and get properties'''

    def __init__(self, props_instance):
        '''define the paramters of the material model'''
        self.parameters = props_instance.get_start_parameters()
        self.mapping = props_instance.get_prop_mapping()
        self.nprops = len(self.parameters)
        self.loaded_data = props_instance

    def get_sorted_props(self):
        '''returns the properties in the order provided in the input props file'''
        return self.loaded_data.get_start_sortet_parameters()

    def get_keys(self):
        '''returns all keys in parameter dictionary'''
        return self.parameters.keys()

    def reset(self, start_props):
        '''quick way to get back to original constants'''
        self.__init__(start_props)

    def getv(self, key):
        '''gets the value of a certain entry'''
        return self.parameters[key]

    def getv_only(self, key):
        '''gets the value of a certain entry'''
        return self.parameters[key][0]

    def setv(self, key, value, ty=float, mini=0, maxi=1, nr=0):
        '''sets a sertian value'''
        self.parameters[key] = (value, ty, mini, maxi, nr)

    def setv_only(self, key, value):
        '''sets a sertian value'''
        val, ty, mini, maxi, nr = self.getv(key)
        self.parameters[key] = (value, ty, mini, maxi, nr)

    def get_umat_input(self):
        '''just returns the sorted list of parameters for the umat call'''
        props = []
        for idx in range(self.nprops):
            props.append(float(self.parameters[self.mapping[idx]][0]))
        return props


class Umat_wrapper():
    '''Class that handles the interfacing with the umat code. Is capable of calling
       fortran code and retrieving stress results from the umat calculaton. It 
       requieres an instance of a material model to be initialized.'''

    def __init__(self, matmod):
        '''init function, creates the instance'''
        self.load_model(matmod)

    def get_props_dic(self):
        '''return the dictionary of properties from the underlying material model'''
        return self.m.parameters

    def setv(self, key, value, ty=float, mini=0, maxi=1, nr=0):
        '''set a value in the underlying material model'''
        self.m.setv(key, value, ty, mini, maxi, nr)

    def getv(self, key):
        '''get a value in the underlying material model'''
        return self.m.getv(key)

    def get_strstres(self):
        '''get the strains and stresses from the loaded experiment'''
        return self.expstrains, self.expstresses

    def plot(self):
        '''after loading plot the data'''
        plt.xlim(1, 1.75)
        plt.ylim(0, 1100)
        for i in range(2):
            plt.plot(self.expstrains[:, i], self.expstresses[:,
                     i], label='Experiment '+str(i), linewidth=3)

    def load_model(self, matmod):
        '''loads in an instance of the material model parameters'''
        self.m = matmod

    def get_model(self):
        '''retrieves material model'''
        return self.m

    def call_umat(self, dfgrd1, props):
        '''calls the umat with the Gcurrent values in the material model'''
        stress = np.array([0. for i in range(6)])
        umat.myumat(stress, props, dfgrd1, len(props))
        return stress

    def load_inputs(self, inp):
        '''loads an input file of strains and stresses'''
        a = np.loadtxt(inp, float, skiprows=1)
        self.expstrains = a[:, :6]
        self.expstresses = a[:, 6:]
        self.dfgrds = [self.convert(strain) for strain in self.expstrains]

    def exp_inputs(self, strains, stresses):
        '''throw in strains and stresses directly'''
        self.expstrains = np.array(strains)
        self.expstresses = np.array(stresses)
        self.dfgrds = [self.convert(strain) for strain in self.expstrains]

    def test_inputs(self, strains):
        '''throw in dformation gradients for a certain test'''
        self.dfgrds = [self.convert(strain) for strain in strains]
        # create experiment stresses and strains based on the default params
        props = self.m.get_umat_input()
        calcstresses = np.array(
            [pk2_voigt(mat6n33a(self.call_umat(df, props)), df) for df in self.dfgrds])
        self.expstresses = np.array(calcstresses)
        self.expstrains = np.array(strains)

    def helper(self, p, *strains):
        '''helper function that solve a system of equations'''
        x11, x12, x13, x22, x23, x33 = p
        strains = strains[0]
        return (x11**2 + x12**2 + x13**2 - strains[0],
                x12**2 + x22**2 + x23**2 - strains[1],
                x13**2 + x23**2 + x33**2 - strains[2],
                x11*x12 + x12*x22 + x13*x23 - strains[3],
                x11*x13 + x12*x23 + x13*x33 - strains[4],
                x12*x13 + x22*x23 + x23*x33 - strains[5]
                )

    def convert(self, strains):
        '''converts strains into deformgradient'''
        # assume that there is no rotation, then F = U, U^2 = C, and C is
        # only diagonal
        x11, x12, x13, x22, x23, x33 =  \
            fsolve(self.helper, (1, 0, 0, 1, 0, 1), args=strains)
        dfgrad = np.zeros((3, 3))
        dfgrad[0][0] = x11
        dfgrad[0][1] = x12
        dfgrad[0][2] = x13
        dfgrad[1][0] = x12
        dfgrad[1][1] = x22
        dfgrad[1][2] = x23
        dfgrad[2][0] = x13
        dfgrad[2][1] = x23
        dfgrad[2][2] = x33
        return dfgrad

    def cost_calculations(self, get_plots=False, axis=[]):
        '''will calculate the cost for a paramter set with a given input,
           if axis = (0) gets cost of 0 axis, can also be (0,1,2)'''
        props = self.m.get_umat_input()
        # the cost is now the difference between the calculated and the
        # given stresses
        calcstresses = np.array(
            [pk2_voigt(mat6n33a(self.call_umat(df, props)), df) for df in self.dfgrds])
        # here I could specify a certain axis to get only the desired cost function of that axis
        if len(axis) == 0:
            print('WARNING! NO COMPONENT SPECIFIED!!!')
            cost = 0
            if get_plots:
                return ([0], [0], 0)

        else:
            cost = 0
            for i in axis:
                cost += np.sum((self.expstresses - calcstresses)**2, 0)[i]
            if get_plots:
                data = []
                for ax in axis:
                    cost = np.sum((self.expstresses - calcstresses)**2, 0)[ax]
                    data.append(
                        (self.expstrains[:, ax], calcstresses[:, ax], cost))
                return data

        return cost


class Load_props():
    '''Class that makes it more convenient to load in a properties file and to handle
       the contents, requieres a path to a filename'''

    def __init__(self, fname):
        '''Initializer'''
        self.data = np.loadtxt(fname, str, skiprows=1)
        self.set_prop_mapping()
        self.nprops = len(self.data)
        self.set_start_parameters()

    def set_prop_mapping(self):
        '''props need to be mapped to entry in props_list for umat call'''
        self.props_map = {}
        for idx in range(len(self.data)):
            self.props_map[idx] = self.data[idx][0]

    def get_prop_mapping(self):
        '''gets the mapping between propertie entry and argument index for umat call'''
        return self.props_map

    def set_start_parameters(self):
        '''retrieves the start values for each of the parameters, order as input'''
        self.prop_dic = {}
        for prop, startval, dtype, mini, maxi in self.data:
            if dtype == 'int':  # for int
                self.prop_dic[prop] = (int(float(startval)), int, int(
                    float(mini)), int(float(maxi)), 0)
            else:
                self.prop_dic[prop] = (
                    float(startval), float, float(mini), float(maxi), 0)

    def get_start_sortet_parameters(self):
        '''returns the start parameters in the order that the umat want'''
        svals = []
        for idx in range(self.nprops):
            prop = self.props_map[idx]
            svals.append((prop, self.prop_dic[prop][0]))
        return svals

    def get_start_parameters(self):
        '''returns the dicitonary of properties'''
        return self.prop_dic


class Convert_experiments():
    '''a class that can deal with biaxial experiments, and writes them 
        into usefull data outputs, for biaxial experiments'''

    def __init__(self):
        pass

    def load_exp(self, filen):
        '''loads in a file'''
        self.data = np.loadtxt(filen, skiprows=1, dtype=float)

    def make_full_strain_stresses(self):
        'calculate the deforamtion fradients from the loaded strains and stresses'
        strains = []
        stresses = []
        dfgrd = np.zeros((3, 3))
        for line in self.data:
            dfgrd[0, 0] = math.sqrt(line[0])
            dfgrd[1, 1] = math.sqrt(line[1])
            dfgrd[2, 2] = 1/math.sqrt(line[0]*line[1])
            c = cauchy_green_strain(dfgrd)
            strains.append([c[0, 0], c[1, 1], c[2, 2],
                           c[0, 1], c[0, 2], c[1, 2]])
            stresses.append([line[2], line[3], 0, 0, 0, 0])
        return np.array(strains), np.array(stresses)


def plot(x, y, idx):
    '''convenience function to show results of an experiment'''
    xs = x[:, idx]
    ys = y[:, idx]
    plt.xlabel('Cauchy Green Strain')
    plt.ylabel('2Piola Kirchhoff Stress')
    plt.title('Curve Fitting for Reese et al Experiments')
    plt.plot(xs, ys, label='umat Calculation '+str(idx), linewidth=3)
