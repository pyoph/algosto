library("paws")

Sys.setenv("AWS_ACCESS_KEY_ID" = "SIOL69ZHL40MFHZ7E219",
           "AWS_SECRET_ACCESS_KEY" = "WhnlPVOpmw8VyrVR6o5ESXUjHL5+yBcdmDOywZ6d",
           "AWS_DEFAULT_REGION" = "us-east-1",
           "AWS_SESSION_TOKEN" = "eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9.eyJhY2Nlc3NLZXkiOiJTSU9MNjlaSEw0ME1GSFo3RTIxOSIsImFjciI6IjAiLCJhbGxvd2VkLW9yaWdpbnMiOlsiKiJdLCJhdWQiOlsibWluaW8iLCJhY2NvdW50Il0sImF1dGhfdGltZSI6MTc1NDkxNjEzNywiYXpwIjoib255eGlhLW1pbmlvIiwiZW1haWwiOiJwYXVsLmd1aWxsb3RAZW5zYWUuZnIiLCJlbWFpbF92ZXJpZmllZCI6dHJ1ZSwiZXhwIjoxNzU2MTI1NzM4LCJmYW1pbHlfbmFtZSI6IkdVSUxMT1QiLCJnaXZlbl9uYW1lIjoiUGF1bCIsImdyb3VwcyI6WyJkYXRhbGFiLWd1aWxsb3QtcGF1bCJdLCJpYXQiOjE3NTQ5MTYxMzgsImlzcyI6Imh0dHBzOi8vYXV0aC5ncm91cGUtZ2VuZXMuZnIvcmVhbG1zL2dlbmVzIiwianRpIjoiOTgyM2JlZTEtYTZlMi00YzliLWE2ZWMtYmI0Yjk5ZjliMjllIiwibmFtZSI6IlBhdWwgR1VJTExPVCIsInBvbGljeSI6InN0c29ubHkiLCJwcmVmZXJyZWRfdXNlcm5hbWUiOiJwZ3VpbGxvdC1lbnNhZSIsInJlYWxtX2FjY2VzcyI6eyJyb2xlcyI6WyJvZmZsaW5lX2FjY2VzcyIsImRlZmF1bHQtcm9sZXMtZ2VuZXMiLCJ1bWFfYXV0aG9yaXphdGlvbiJdfSwicmVzb3VyY2VfYWNjZXNzIjp7ImFjY291bnQiOnsicm9sZXMiOlsibWFuYWdlLWFjY291bnQiLCJtYW5hZ2UtYWNjb3VudC1saW5rcyIsInZpZXctcHJvZmlsZSJdfX0sInNjb3BlIjoib3BlbmlkIHByb2ZpbGUgZW1haWwiLCJzaWQiOiIwNDlmMTQwMi0zMTFjLTQwZTQtOWE2My0zMjgwZTY4ZDgyODgiLCJzdWIiOiI3OTVlMjFiNC1iZmI1LTRmZjAtYjViMS1lNDNlYWE3N2NkMTAiLCJ0eXAiOiJCZWFyZXIifQ.RCRy18UAljjj3rQ6U5b9Ko0mxpxWxFLhZbiMco4PNxe-Oko9yqT1zl9SX0cuYFrxFLbOZ9dnWiCch7VQv0Z38w",
           "AWS_S3_ENDPOINT"= "minio-simple.lab.groupe-genes.fr")

minio <- paws::s3(config = list(
  credentials = list(
    creds = list(
      access_key_id = Sys.getenv("AWS_ACCESS_KEY_ID"),
      secret_access_key = Sys.getenv("AWS_SECRET_ACCESS_KEY"),
      session_token = Sys.getenv("AWS_SESSION_TOKEN")
    )),
  endpoint = paste0("https://", Sys.getenv("AWS_S3_ENDPOINT")),
  region = Sys.getenv("AWS_DEFAULT_REGION")))

minio$list_buckets()

bucket_name <- "projet-datalab-guillot-paul"
key <- "FitSimKrho/FitSim/FitParms-d10-n10000-k0.858565436437754-l1-rho0.85-r25-sim75.RData"

# Télécharger le fichier
obj <- minio$get_object(Bucket = bucket_name, Key = key)

# Sauvegarder dans un fichier temporaire
tmp <- tempfile(fileext = ".RData")
writeBin(obj$Body, tmp)

# Charger dans l'environnement courant
load(tmp, envir = .GlobalEnv)

message("Fichier chargé : ", key)