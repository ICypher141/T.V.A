import spacy
from astroquery.simbad import Simbad
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
from astropy import units as u
from datetime import datetime, timezone
from math import sin, cos, asin, radians, degrees

# Load the spaCy language model
nlp = spacy.load("custom_celestial_ner")

# Function to fetch celestial object coordinates
def fetch_celestial_object_coordinates(object_name):
    try:
        result_table = Simbad.query_object(object_name)
        if result_table:
            ra = result_table['RA'][0]  # Right Ascension
            dec = result_table['DEC'][0]  # Declination
            return {"RA": ra, "Dec": dec}
        else:
            return None
    except Exception as e:
        return f"Error: {e}"

# Function to parse commands
def parse_command(command):
    doc = nlp(command)
    intent = None
    object_name = None

    # Extract intent based on keywords
    for token in doc:
        if token.lemma_ in ["find", "point", "show", "give"]:
            intent = token.lemma_

    # Extract celestial object names (proper nouns)
    for ent in doc.ents:
        if ent.label_ == "CELESTIAL_OBJECT":
            object_name = ent.text

    return intent, object_name

# Function to convert RA and Dec to degrees
def convert_ra_dec_to_degrees(ra, dec):
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
    return coord.ra.deg, coord.dec.deg

# Function to calculate current zenith
def calculate_zenith(location):
    observing_time = Time.now()
    altaz_frame = AltAz(location=location, obstime=observing_time)
    zenith = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=altaz_frame).transform_to("icrs")
    return {"RA": zenith.ra.to_string(unit=u.hour), "Dec": zenith.dec.to_string(unit=u.deg)}

# Function to calculate altitude and azimuth
def calculate_alt_az(latitude, longitude, ra_deg, dec_deg):
    current_time = datetime.now(timezone.utc)

    # Local Sidereal Time (LST)
    JD = (current_time - datetime(2000, 1, 1, 12, tzinfo=timezone.utc)).total_seconds() / 86400 + 2451545.0
    T = (JD - 2451545.0) / 36525.0
    GMST = 280.46061837 + 360.98564736629 * (JD - 2451545.0) + T**2 * (0.000387933 - T / 38710000.0)
    GMST %= 360.0
    lst = (GMST + longitude) % 360.0

    # Hour Angle
    hour_angle = (lst - ra_deg) % 360.0
    hour_angle = hour_angle if hour_angle <= 180 else hour_angle - 360

    # Altitude and Azimuth
    lat_rad = radians(latitude)
    dec_rad = radians(dec_deg)
    ha_rad = radians(hour_angle)

    altitude = degrees(asin(sin(lat_rad) * sin(dec_rad) + cos(lat_rad) * cos(dec_rad) * cos(ha_rad)))

    az_sin = -cos(lat_rad) * sin(ha_rad) / cos(radians(altitude))
    azimuth = degrees(asin(az_sin))

    if hour_angle > 0:
        azimuth = 180 - azimuth
    else:
        azimuth = 180 + azimuth

    return {"Altitude": altitude, "Azimuth": azimuth}
