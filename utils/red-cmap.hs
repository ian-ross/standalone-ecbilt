module Main where

default (Int, Double)

main :: IO ()
main = do
  let delta1 = (255.0 - 10.0) / 126.0 :: Double
      delta2 = (255.0 - 128.0) / 126.0 :: Double
      s1 i = "255 " ++ show i ++ " " ++ show i
      i1 = map (fromIntegral . round) [255, 255-delta1 .. 10]
      pt1 = map s1 i1
      s2 i = show i ++ " 0 0"
      i2 = map (fromIntegral . round) [255, 255-delta2 .. 128]
      pt2 = map s2 i2
  putStrLn $ "ncolors=" ++ show (length pt1 + length pt2)
  mapM_ putStrLn $ pt1 ++ pt2
