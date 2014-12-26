classdef Molecule < handle
    
    properties (SetAccess = private)
        
        cartesian;
        charge = 0;
        multiplicity = 1;
        
    end
        
    methods
        
        function obj = Molecule(cartesian, charge, multiplicity)
            if(nargin > 1)
                obj.charge = charge;
            end
            if(nargin > 2)
                obj.multiplicity = multiplicity;
            end
            obj.cartesian = cartesian;
        end
        
        function molString = MoleculeString(obj)
            numAtoms = size(obj.cartesian,1);
            symbols = repmat(' ', numAtoms, 5);
            for iatom = 1:numAtoms
                atomSymbol = obj.PeriodicTable{obj.cartesian(iatom, 1)};
                symbols(iatom, 1:length(atomSymbol)) = atomSymbol;
            end
            formSpec = '%.10f';
            spaces = repmat('  ', numAtoms, 1);
            molString = [symbols, num2str(obj.cartesian(:, 2), formSpec), ...
                spaces, num2str(obj.cartesian(:, 3), formSpec), ...
                spaces, num2str(obj.cartesian(:, 4), formSpec), ...
                repmat(char(10), numAtoms, 1)];
            molString = reshape(molString', 1, []);
        end
        
        function numAtoms = NumAtoms(obj)
            numAtoms = size(obj.cartesian, 1);
        end
        
        function numElectrons = NumElectrons(obj)
            numElectrons = sum(obj.cartesian(:, 1)) - obj.charge;
        end
        
        function distMat = DistanceMatrix(obj)
            distMat = dist(obj.cartesian(:, 2:end)');
        end
        
    end
    
    methods (Static)
        
        cartesian = ZMatrixToCartesian(zMatrix);
        zmatrix = ZMatStrToZMatrix(zMatStr);
        
    end
    
    properties (Constant, Access = private)
        
        PeriodicTable = Molecule.InitializePeriodicTable()
        Bohr2Angstrom = 0.529177249;
        
    end
    
    methods (Access = private, Static)
        
        function periodicTable = InitializePeriodicTable()
            periodicTable = { ...
                'H',                                      'He', ...
                'Li', 'Be',  'B',  'C',  'N',  'O',  'F', 'Ne', ...
                'Na', 'Mg', 'Al', 'Si',  'P',  'S', 'Cl', 'Ar', ...
                'K',  'Ca'};
        end
        
    end
    
end