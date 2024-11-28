import streamlit as st
from model import fetch_celestial_object_coordinates, parse_command, convert_ra_dec_to_degrees, calculate_zenith, calculate_alt_az
from astropy.coordinates import EarthLocation
import speech_recognition as sr
import time
from astropy import units as u
import sounddevice as sd

# Observer's location
latitude = 18.5204
longitude = 73.8567
location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg, height=560 * u.m)

# Global storage for celestial object coordinates and results
st.session_state.celestial_coordinates = st.session_state.get("celestial_coordinates", {})
st.session_state.alt_az = st.session_state.get("alt_az", None)
st.session_state.zenith_radec = st.session_state.get("zenith_radec", None)

# App Title
st.title("Celestial Object Tracker")

# Step 1: Capture Voice Command
st.header("Capture Command and Process Automatically")

if st.button("Record and Process"):
    recognizer = sr.Recognizer()
    with sr.Microphone() as source:
        st.write("Adjusting for background noise...")
        recognizer.adjust_for_ambient_noise(source)
        st.write("Listening for your command...")
        try:
            # Capture and process the audio
            audio = recognizer.listen(source, timeout=5, phrase_time_limit=5)
            st.write("Processing your command...")
            command = recognizer.recognize_google(audio)
            st.success(f"Recognized command: {command}")

            # Parse the command
            intent, object_name = parse_command(command)
            if intent and object_name:
                st.write(f"Intent: {intent}, Object Name: {object_name}")

                # Fetch celestial coordinates
                coordinates = fetch_celestial_object_coordinates(object_name)
                if coordinates:
                    # Convert to degrees
                    ra_deg, dec_deg = convert_ra_dec_to_degrees(coordinates["RA"], coordinates["Dec"])
                    st.session_state.celestial_coordinates = {
                        "RA_deg": ra_deg,
                        "Dec_deg": dec_deg,
                        "RA": coordinates["RA"],
                        "Dec": coordinates["Dec"],
                        "Object Name": object_name,
                    }
                    st.success(
                        f"Coordinates for {object_name}: RA (degrees): {ra_deg}, Dec (degrees): {dec_deg}"
                    )

                    # Calculate Alt/Az
                    st.session_state.alt_az = calculate_alt_az(
                        latitude, longitude, ra_deg, dec_deg
                    )
                else:
                    st.error("Object not found in Simbad.")
            else:
                st.error("Could not parse command.")
        except sr.WaitTimeoutError:
            st.error("Timeout: No speech detected.")
        except sr.UnknownValueError:
            st.error("Sorry, I could not understand the audio.")
        except sr.RequestError as e:
            st.error(f"Request error: {e}")

    # Step 2: Calculate current zenith RA/Dec
    st.session_state.zenith_radec = calculate_zenith(location)

# Display results dynamically
st.header("Results")

if st.session_state.zenith_radec:
    zenith = st.session_state.zenith_radec
    st.subheader("Current Zenith RA/Dec")
    st.write(f"RA: {zenith['RA']}")
    st.write(f"Dec: {zenith['Dec']}")

if st.session_state.celestial_coordinates:
    celestial = st.session_state.celestial_coordinates
    st.subheader(f"Celestial Object: {celestial['Object Name']}")
    st.write(f"RA (degrees): {celestial['RA_deg']}")
    st.write(f"Dec (degrees): {celestial['Dec_deg']}")

if st.session_state.alt_az:
    altaz = st.session_state.alt_az
    st.subheader("Calculated Altitude and Azimuth")
    st.write(f"Altitude: {altaz['Altitude']:.2f}°")
    st.write(f"Azimuth: {altaz['Azimuth']:.2f}°")
