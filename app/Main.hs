{-# LANGUAGE OverloadedStrings #-}

module Main where

import Data.Aeson
import Data.ByteString.Lazy (ByteString)
import qualified Data.ByteString.Lazy.Char8 as BL
import System.IO (hSetBuffering, stdout, stdin, BufferMode(..))

import LimCalc.Core.Expr
import LimCalc.Core.ExprJSON ()
import LimCalc.Core.Simplify
import LimCalc.Pretty
import LimCalc.Differentiation.Calculus
import LimCalc.Differentiation.Limit
import LimCalc.Integration.Risch

-- | A JSON command sent from the Python client.
data Command = Command
  { cmdOp   :: String
  , cmdExpr :: Maybe Expr
  , cmdVar  :: Maybe String
  , cmdX0   :: Maybe Double
  , cmdVars :: Maybe [String]
  } deriving (Show)

instance FromJSON Command where
  parseJSON = withObject "Command" $ \o -> Command
    <$> o .:  "op"
    <*> o .:? "expr"
    <*> o .:? "var"
    <*> o .:? "x0"
    <*> o .:? "vars"

-- | A JSON response sent back to the Python client.
data Response
  = ROk Value
  | RError String

instance ToJSON Response where
  toJSON (ROk v)    = object ["ok" .= True,  "result" .= v]
  toJSON (RError e) = object ["ok" .= False, "error"  .= e]

-- | Dispatch a command to the appropriate limcalc operation.
dispatch :: Command -> Response
dispatch cmd = case cmdOp cmd of

  "diff" ->
    withExpr cmd $ \f ->
    withVar  cmd $ \v ->
    case diff f v of
      Left err -> RError (show err)
      Right e  -> ROk (toJSON (simplify e))

  "partial_diff" ->
    withExpr cmd $ \f ->
    withVar  cmd $ \v ->
    case partialDiff f v of
      Left err -> RError (show err)
      Right e  -> ROk (toJSON (simplify e))

  "integrate" ->
    withExpr cmd $ \f ->
    withVar  cmd $ \v ->
    case rischIntegrate f v of
      Elementary e   -> ROk (toJSON e)
      NonElementary  -> RError "NonElementary"
      NotImplemented s -> RError ("NotImplemented: " ++ s)
      RischError s   -> RError ("RischError: " ++ s)

  "limit" ->
    withExpr cmd $ \f ->
    withVar  cmd $ \v ->
    withX0   cmd $ \x0 ->
    case limit f v x0 of
      Exists val     -> ROk (toJSON val)
      Pole r         -> RError ("Pole at order " ++ show r)
      DoesNotExist s -> RError ("DoesNotExist: " ++ s)
      LimitError e   -> RError (show e)

  "simplify" ->
    withExpr cmd $ \f ->
    ROk (toJSON (simplify f))

  "pretty" ->
    withExpr cmd $ \f ->
    ROk (toJSON (prettyExpr f))

  "gradient" ->
    withExpr cmd $ \f ->
    withVars cmd $ \vs ->
    case gradient f vs of
      Left err -> RError (show err)
      Right es -> ROk (toJSON (map simplify es))

  op -> RError ("Unknown operation: " ++ op)

-- Helpers for extracting required fields

withExpr :: Command -> (Expr -> Response) -> Response
withExpr cmd f = case cmdExpr cmd of
  Nothing -> RError "Missing 'expr' field"
  Just e  -> f e

withVar :: Command -> (String -> Response) -> Response
withVar cmd f = case cmdVar cmd of
  Nothing -> RError "Missing 'var' field"
  Just v  -> f v

withVars :: Command -> ([String] -> Response) -> Response
withVars cmd f = case cmdVars cmd of
  Nothing -> RError "Missing 'vars' field"
  Just vs -> f vs

withX0 :: Command -> (Double -> Response) -> Response
withX0 cmd f = case cmdX0 cmd of
  Nothing -> RError "Missing 'x0' field"
  Just x  -> f x

-- | Main loop: read one JSON command per line, write one JSON response per line.
main :: IO ()
main = do
  hSetBuffering stdout LineBuffering
  hSetBuffering stdin  LineBuffering
  loop

loop :: IO ()
loop = do
  line <- BL.pack <$> getLine
  case eitherDecode line of
    Left err  -> BL.putStrLn (encode (RError ("Parse error: " ++ err)))
    Right cmd -> BL.putStrLn (encode (dispatch cmd))
  loop
