-- | JSON serialization and deserialization for 'Expr'.
--
-- Provides 'ToJSON' and 'FromJSON' instances for 'Expr', enabling
-- the limcalc engine to communicate with external clients (Python,
-- REST APIs, CLI tools) via JSON.
--
-- = JSON representation
--
-- Each constructor is represented as a JSON object with a @\"tag\"@
-- field identifying the constructor, plus constructor-specific fields:
--
-- @
-- Const 2.0        → {\"tag\": \"Const\", \"value\": 2.0}
-- Var \"x\"          → {\"tag\": \"Var\", \"name\": \"x\"}
-- Pi               → {\"tag\": \"Pi\"}
-- E                → {\"tag\": \"E\"}
-- I                → {\"tag\": \"I\"}
-- Add f g          → {\"tag\": \"Add\", \"left\": ..., \"right\": ...}
-- Sub f g          → {\"tag\": \"Sub\", \"left\": ..., \"right\": ...}
-- Mul f g          → {\"tag\": \"Mul\", \"left\": ..., \"right\": ...}
-- Div f g          → {\"tag\": \"Div\", \"left\": ..., \"right\": ...}
-- Pow f g          → {\"tag\": \"Pow\", \"base\": ..., \"exp\": ...}
-- Neg f            → {\"tag\": \"Neg\", \"arg\": ...}
-- Abs f            → {\"tag\": \"Abs\", \"arg\": ...}
-- Exp f            → {\"tag\": \"Exp\", \"arg\": ...}
-- Log f            → {\"tag\": \"Log\", \"arg\": ...}
-- Sin f            → {\"tag\": \"Sin\", \"arg\": ...}
-- Cos f            → {\"tag\": \"Cos\", \"arg\": ...}
-- Arcsin f         → {\"tag\": \"Arcsin\", \"arg\": ...}
-- Arccos f         → {\"tag\": \"Arccos\", \"arg\": ...}
-- Arctan f         → {\"tag\": \"Arctan\", \"arg\": ...}
-- Erf f            → {\"tag\": \"Erf\", \"arg\": ...}
-- Li f             → {\"tag\": \"Li\", \"arg\": ...}
-- Si f             → {\"tag\": \"Si\", \"arg\": ...}
-- Ci f             → {\"tag\": \"Ci\", \"arg\": ...}
-- Ei f             → {\"tag\": \"Ei\", \"arg\": ...}
-- @

{-# LANGUAGE OverloadedStrings #-}

module LimCalc.Core.ExprJSON
  ( -- * Re-exports
    module LimCalc.Core.Expr
  ) where

import Data.Aeson
import Data.Aeson.Types (Parser)
import LimCalc.Core.Expr

instance ToJSON Expr where
  toJSON (Const d)    = object ["tag" .= ("Const" :: String), "value" .= d]
  toJSON Pi           = object ["tag" .= ("Pi" :: String)]
  toJSON E            = object ["tag" .= ("E" :: String)]
  toJSON I            = object ["tag" .= ("I" :: String)]
  toJSON (Var x)      = object ["tag" .= ("Var" :: String), "name" .= x]
  toJSON (Add f g)    = object ["tag" .= ("Add" :: String), "left" .= f, "right" .= g]
  toJSON (Sub f g)    = object ["tag" .= ("Sub" :: String), "left" .= f, "right" .= g]
  toJSON (Mul f g)    = object ["tag" .= ("Mul" :: String), "left" .= f, "right" .= g]
  toJSON (Div f g)    = object ["tag" .= ("Div" :: String), "left" .= f, "right" .= g]
  toJSON (Pow f g)    = object ["tag" .= ("Pow" :: String), "base" .= f, "exp" .= g]
  toJSON (Neg f)      = object ["tag" .= ("Neg" :: String), "arg" .= f]
  toJSON (Abs f)      = object ["tag" .= ("Abs" :: String), "arg" .= f]
  toJSON (Exp f)      = object ["tag" .= ("Exp" :: String), "arg" .= f]
  toJSON (Log f)      = object ["tag" .= ("Log" :: String), "arg" .= f]
  toJSON (Sin f)      = object ["tag" .= ("Sin" :: String), "arg" .= f]
  toJSON (Cos f)      = object ["tag" .= ("Cos" :: String), "arg" .= f]
  toJSON (Arcsin f)   = object ["tag" .= ("Arcsin" :: String), "arg" .= f]
  toJSON (Arccos f)   = object ["tag" .= ("Arccos" :: String), "arg" .= f]
  toJSON (Arctan f)   = object ["tag" .= ("Arctan" :: String), "arg" .= f]
  toJSON (Erf f)      = object ["tag" .= ("Erf" :: String), "arg" .= f]
  toJSON (Li f)       = object ["tag" .= ("Li" :: String), "arg" .= f]
  toJSON (Si f)       = object ["tag" .= ("Si" :: String), "arg" .= f]
  toJSON (Ci f)       = object ["tag" .= ("Ci" :: String), "arg" .= f]
  toJSON (Ei f)       = object ["tag" .= ("Ei" :: String), "arg" .= f]

instance FromJSON Expr where
  parseJSON = withObject "Expr" $ \o -> do
    tag <- o .: "tag" :: Parser String
    case tag of
      "Const"  -> Const  <$> o .: "value"
      "Pi"     -> pure Pi
      "E"      -> pure E
      "I"      -> pure I
      "Var"    -> Var    <$> o .: "name"
      "Add"    -> Add    <$> o .: "left"  <*> o .: "right"
      "Sub"    -> Sub    <$> o .: "left"  <*> o .: "right"
      "Mul"    -> Mul    <$> o .: "left"  <*> o .: "right"
      "Div"    -> Div    <$> o .: "left"  <*> o .: "right"
      "Pow"    -> Pow    <$> o .: "base"  <*> o .: "exp"
      "Neg"    -> Neg    <$> o .: "arg"
      "Abs"    -> Abs    <$> o .: "arg"
      "Exp"    -> Exp    <$> o .: "arg"
      "Log"    -> Log    <$> o .: "arg"
      "Sin"    -> Sin    <$> o .: "arg"
      "Cos"    -> Cos    <$> o .: "arg"
      "Arcsin" -> Arcsin <$> o .: "arg"
      "Arccos" -> Arccos <$> o .: "arg"
      "Arctan" -> Arctan <$> o .: "arg"
      "Erf"    -> Erf    <$> o .: "arg"
      "Li"     -> Li     <$> o .: "arg"
      "Si"     -> Si     <$> o .: "arg"
      "Ci"     -> Ci     <$> o .: "arg"
      "Ei"     -> Ei     <$> o .: "arg"
      _        -> fail $ "Unknown Expr tag: " ++ tag
