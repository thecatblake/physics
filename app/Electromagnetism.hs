module Electromagnetism where
    import Common
    import Linear.Vector
    import Linear.Metric

    coulombConstant :: Double
    coulombConstant = 8.99e9

    electronCoulomb :: Double
    electronCoulomb = 1.602e-19

    electronMass :: Double
    electronMass = 9.1093837e-31

    

    data PointCharge = PointCharge RVec3 Double

    data VoltageUnit =
        Volt Double |
        EVolt Double
        deriving (Eq, Show)
    
    voltUnitToVolt :: VoltageUnit -> Double
    voltUnitToVolt (Volt x) = x
    voltUnitToVolt (EVolt x) = electronCoulomb * x

    coulombForce' :: Double -> Double -> Double -> Double
    coulombForce' q q' r = coulombConstant * q * q' / (r * r)

    coulombForce :: PointCharge -> PointCharge -> RVec3
    coulombForce (PointCharge r q) (PointCharge r' q') = (coulombForce' q q' l2) *^ r
                                                            where
                                                                l2 = norm (r - r')