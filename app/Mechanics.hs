module Mechanics where
    import Common 
    import Linear.Vector
    import Linear.V3
    import Linear.Metric

    data MassUnit = 
        Gram Double |
        KiloGram Double |
        Pounds Double 
        deriving (Eq, Ord, Show)

    data TimeUnit = 
        MiliSeconds Double |
        Seconds Double
        deriving (Eq, Ord, Show)

    data PointMass = PointMass Double RVec3 RVec3 deriving (Eq, Show)

    gravityConstant = 9.7803
    
    massToKilo :: MassUnit -> Double
    massToKilo (Gram x) = x / 1000
    massToKilo (KiloGram x) = x
    massToKilo (Pounds x) = 0.453592 * x

    massToPounds :: MassUnit -> Double
    massToPounds (Gram x) = 0.00220462 * x
    massToPounds (KiloGram x) = 2.20462 * x
    massToPounds (Pounds x) = x

    timeToSeconds :: TimeUnit -> Double
    timeToSeconds (MiliSeconds x) = x / 1000
    timeToSeconds (Seconds x) = x

    kineticEnergy :: PointMass -> Double
    kineticEnergy (PointMass m r v) = m * v' * v' / 2
                                        where
                                            v' = norm v

    elasticCollision :: PointMass -> PointMass -> (PointMass, PointMass)
    elasticCollision (PointMass m1 r1 v1) (PointMass m2 r2 v2) = (PointMass m1 r1 ((m'/m) *^ v1 + (2*m2 / m) *^ v2), PointMass m2 r2 (((-m')/m) *^ v2 + (2*m1 / m) *^ v1))
                                                                where 
                                                                    m = m1 + m2
                                                                    m' = m1 - m2

    momentum :: PointMass -> RVec3
    momentum (PointMass m _ v) = m *^ v
    
    centerOfMass :: [PointMass] -> PointMass
    centerOfMass = foldr (\(PointMass m1 r1 v1) (PointMass m2 r2 v2) -> let m = m1 + m2 in PointMass m ((1/m) *^ (m1 *^ r1 + m2 *^ r2)) ((1/m) *^ (m1 *^ v1 + m2 *^ v2))) (PointMass 0 zeroV zeroV)

    freeFallDist :: Double -> Double
    freeFallDist t = gravityConstant * t * t / 2

    freeFallTime :: Double -> Double
    freeFallTime d = sqrt (2 * d / gravityConstant)

    kineticToGravityPotential :: PointMass -> PointMass
    kineticToGravityPotential (PointMass m r v) = PointMass m ((V3 0 ((norm v)^2 / (2 * gravityConstant)) 0) + r) zeroV