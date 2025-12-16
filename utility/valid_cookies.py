class ValidCookies:
    def middleware(self, site_cookies: list)
        def validate(request):
            cookies = getattr(request, 'cookies', {})
            for element in site_cookies: 
                if element not in cookies or not cookies[name]:
                    return {"Alert": f"Missing cookie: {element}"}, 401 
            return None
        return validate

