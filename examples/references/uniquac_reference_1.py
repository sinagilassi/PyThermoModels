# ====================================
# CUSTOM REFERENCES
# ====================================
REFERENCE_CONTENT = """
REFERENCES:
    CUSTOM-REFERENCE-1:
      DATABOOK-ID: 1
      TABLES:
        UNIQUAC-binary-parameters:
            TABLE-ID: 1
            DESCRIPTION:
                Directional UNIQUAC binary interaction-energy parameters transcribed from
                the supplied table. Matrix direction is row component i to column component j.
                Diagonal dU values are set to zero. The source image does not identify the
                energy unit; replace SOURCE-UNIT-UNRESOLVED with the confirmed unit before use.
            MATRIX-SYMBOL:
                - UNIQUAC interaction-energy difference: dU
            STRUCTURE:
                COLUMNS: [No., Mixture, Name, Formula, State, dU_i_1, dU_i_2]
                SYMBOL: [None, None, None, None, None, dU_i_1, dU_i_2]
                UNIT: [None, None, None, None, None, SOURCE-UNIT-UNRESOLVED, SOURCE-UNIT-UNRESOLVED]
            VALUES:
                - [1, methanol|ethanol, methanol, CH3OH, l, 0, -1185.15]
                - [2, methanol|ethanol, ethanol, C2H5OH, l, 1670.482, 0]
                - [1, methanol|dimethyl-carbonate, methanol, CH3OH, l, 0, -191.038]
                - [2, methanol|dimethyl-carbonate, dimethyl-carbonate, C3H6O3, l, 2620.497, 0]
                - [1, methanol|ethyl-methyl-carbonate, methanol, CH3OH, l, 0, 47.87923]
                - [2, methanol|ethyl-methyl-carbonate, ethyl-methyl-carbonate, C4H8O3, l, 2271.729, 0]
                - [1, methanol|diethyl-carbonate, methanol, CH3OH, l, 0, -462.978]
                - [2, methanol|diethyl-carbonate, diethyl-carbonate, C5H10O3, l, 3108.645, 0]
                - [1, ethanol|dimethyl-carbonate, ethanol, C2H5OH, l, 0, -440.509]
                - [2, ethanol|dimethyl-carbonate, dimethyl-carbonate, C3H6O3, l, 2415.906, 0]
                - [1, ethanol|ethyl-methyl-carbonate, ethanol, C2H5OH, l, 0, -156.138]
                - [2, ethanol|ethyl-methyl-carbonate, ethyl-methyl-carbonate, C4H8O3, l, 1726.915, 0]
                - [1, ethanol|diethyl-carbonate, ethanol, C2H5OH, l, 0, 691.9057]
                - [2, ethanol|diethyl-carbonate, diethyl-carbonate, C5H10O3, l, 785.9151, 0]
                - [1, dimethyl-carbonate|ethyl-methyl-carbonate, dimethyl-carbonate, C3H6O3, l, 0, 36.89908]
                - [2, dimethyl-carbonate|ethyl-methyl-carbonate, ethyl-methyl-carbonate, C4H8O3, l, 7.335787, 0]
                - [1, dimethyl-carbonate|diethyl-carbonate, dimethyl-carbonate, C3H6O3, l, 0, 2533.673]
                - [2, dimethyl-carbonate|diethyl-carbonate, diethyl-carbonate, C5H10O3, l, -1604.67, 0]
                - [1, ethyl-methyl-carbonate|diethyl-carbonate, ethyl-methyl-carbonate, C4H8O3, l, 0, -5.19434]
                - [2, ethyl-methyl-carbonate|diethyl-carbonate, diethyl-carbonate, C5H10O3, l, -35.3918, 0]
        general-data:
            TABLE-ID: 2
            DESCRIPTION:
                This table provides pure component UNIQUAC volume and surface-area
                parameters.
            DATA: []
            STRUCTURE:
                COLUMNS: [No.,Name,Formula,State,Molecular-Weight,Critical-Temperature,Critical-Pressure,Critical-Molar-Volume,Critical-Compressibility-Factor,Acentric-Factor,Enthalpy-of-Formation,Gibbs-Energy-of-Formation,Volume-Parameter,Surface-Area-Parameter]
                SYMBOL: [None,None,None,None,MW,Tc,Pc,Vc,Zc,AcFa,EnFo,GiEnFo,r,q]
                UNIT: [None,None,None,None,g/mol,K,MPa,m3/kmol,None,None,kJ/mol,kJ/mol,None,None]
                CONVERSION: [None,None,None,None,1,1,1,1,1,1,1,1,1,1]
            VALUES:
                - [1,'methanol','CH3OH','l',32.04,512.5,8.084,0.117,0.222,0.5658,-200.7,-162,1.4311,1.4320]
                - [2,'ethanol','C2H5OH','l',46.068,514,6.137,0.168,0.241,0.6436,-277.70,-174.80,2.1055,1.8920]
                - [3,'butyl-methyl-ether','C5H12O','l',88.15,497.15,3.43,0.329,0.273,0.2662,-283.5,-117.5,4.0672,3.4920]
                - [4,'dimethyl-carbonate','C3H6O3','l',90.078,557.0,4.908865,0.252,0.267,0.346,None,None,3.0613,2.8160]
                - [5,'ethyl-methyl-carbonate','C4H8O3','l',104.10,546.24,3.83878,0.316,0.267,0.4010,None,None,3.7357,3.3560]
                - [6,'diethyl-carbonate','C5H10O3','l',118.13,578.10,3.41986,0.375,0.267,0.4531,None,None,4.4101,3.8960]
"""
