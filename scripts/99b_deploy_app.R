install.packages("rsconnect")

# Authenticate via browser (opens a login/authorize page) instead of the
# old copy-paste token/secret from setAccountInfo()
rsconnect::connectCloudUser()

# Confirm it worked:
rsconnect::accounts()

# Deploy (same call as before, pointing at your app's folder)
rsconnect::deployApp(appDir = "interactive_app/")